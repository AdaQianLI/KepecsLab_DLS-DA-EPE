#include <Arduino.h>
#include <FastLED.h>
#include "src/MotorShield/MotorShield.h"

/*
  Joystick_ForceField_v3
  ----------------------
  Cleaned release-oriented refactor of the original joystick force-field control
  sketch used with an Arduino Mega 2560.

  Purpose
  - Sample joystick x/y position at 1 kHz using Timer5 + ADC interrupts.
  - Receive force commands from Bpod over Serial1.
  - Drive Motor A continuously during Bpod-defined force states.
  - Drive Motor B for the slide-close / recenter sequence.
  - Return brief TTL pulses to Bpod for force-on and too-early event detection.
  - Provide a reach-level TTL signal, buzzer feedback, and LED feedback.

  Behavioral note
  - Motor force is applied continuously while the joystick is above the
    force-on threshold and FORCE_ENABLE_PIN remains high.
  - TTL_FORCE_ON_PIN returns repeated brief pulses to Bpod while the joystick
    remains above the force-on threshold and FORCE_ONSET_MONITOR_PIN is high.
    This repeated pulsing is intentional so that Bpod does not miss a single
    brief event pulse.
  - TTL_TOO_EARLY_PIN returns repeated brief pulses to Bpod while the
    too-early criterion is satisfied. This repeated pulsing is also intentional.
  - TTL_REACH_PIN is a sustained level signal while the reach criterion is met.
*/


// -----------------------------------------------------------------------------
// Build options
// -----------------------------------------------------------------------------
constexpr bool kEnableSerialDebug = true;   // USB Serial Plotter debug output
constexpr uint32_t kSystemBaudRate = 115200; // Use fast serial if debug is enabled
constexpr uint32_t kBpodBaudRate = 115200;

// -----------------------------------------------------------------------------
// Hardware / task configuration
// -----------------------------------------------------------------------------
constexpr uint8_t kNumLeds = 72;
constexpr uint8_t kLedBrightness = 10;
constexpr uint16_t kAdcPeriodUs = 1000;      // 1 kHz joystick sampling
constexpr uint16_t kTtlPulseMs = 2;          // Pulse width for Bpod return pulses
constexpr uint16_t kMotorCommandPauseUs = 300;

// Joystick baseline and thresholds (ADC units)
volatile int16_t gBaselineX = 354;
volatile int16_t gBaselineY = 368;
constexpr int16_t kXWindowLimit = 60;
constexpr int16_t kForceOffsetThresholdY = 50;  // reach / force-off threshold
constexpr int16_t kForceOnsetThresholdY = 20;   // force-on threshold
constexpr int16_t kTooEarlyThresholdY = 100;    // early movement threshold

// Buzzer shaping thresholds (ADC units)
constexpr int16_t kBuzzerXGate = 50;
constexpr int16_t kBuzzerYGate = 40;
constexpr uint16_t kBaseBuzzerHz = 2000;

// -----------------------------------------------------------------------------
// GPIO mapping
// -----------------------------------------------------------------------------
constexpr uint8_t kLedDataOutPin = 37;   // WS2811 data pin
constexpr uint8_t kLedSeq1TrigPin = 22;  // Bpod BNC OUT 1
constexpr uint8_t kLedSeq2TrigPin = 23;  // Bpod BNC OUT 2
constexpr uint8_t kSlideTrigPin = 2;     // Bpod Wire OUT 1
constexpr uint8_t kBuzzerTrigPin = 24;   // Bpod Wire OUT 2; state gate for position monitoring and buzzer output
constexpr uint8_t kForceOnsetMonitorPin = 25;   // Bpod Wire OUT 3
constexpr uint8_t kForceEnablePin = 26;  // Bpod Wire OUT 4
constexpr uint8_t kBuzzerOutPin = 36;
constexpr uint8_t kTtlReachPin = 30;     // Bpod Wire IN 1
constexpr uint8_t kTtlForceOnPin = 31;   // Bpod Wire IN 2
constexpr uint8_t kTtlTooEarlyPin = 32;    // Bpod Wire IN 3
constexpr uint8_t kPosXPin = A8;
constexpr uint8_t kPosYPin = A9;
constexpr uint8_t kInternalLedPin = 13;

// -----------------------------------------------------------------------------
// Timer / ADC constants
// -----------------------------------------------------------------------------
constexpr uint16_t kTimerMax = 0xFFFF;

// -----------------------------------------------------------------------------
// Global state shared across ISRs and loop()
// -----------------------------------------------------------------------------
volatile int16_t gPosX = 0;
volatile int16_t gPosY = 0;
volatile uint8_t gCurrentAdcChannel = kPosXPin;
volatile bool gAdcSampleReady = false;
volatile bool gSlideResetRequested = false;
volatile bool gBaselineResetRequested = false;

// Bpod command byte that controls motor amplitude
uint8_t gForceFactor = 0;

// Diagnostics / derived joystick variables
int16_t gDiffX = 0;
int16_t gDiffY = 0;

// Peripheral objects
CRGB gLeds[kNumLeds];
MotorShield gMotorShield;

// -----------------------------------------------------------------------------
// Forward declarations
// -----------------------------------------------------------------------------
uint16_t calculateTimerTop(const uint32_t fosc, const uint16_t prescaler, const uint32_t timeTopUs);
void initTimer5(const uint32_t microSeconds);
void initADC();
void handleSlideTrigger();

void readBpodForceCommand();
void processLatestJoystickSample();
void processJoystickSample(const int16_t posX, const int16_t posY);
void updateForceMotor(const int16_t diffX, const int16_t diffY);
void updateTooEarlyPulse(const int16_t diffX, const int16_t diffY);
void updateBuzzerAndReachTtl(const int16_t diffX, const int16_t diffY);
void updateLedSequences();
void recalibrateBaseline();
void runSlideCloseSequence();
void sendPulse(const uint8_t pin, const uint16_t pulseMs);
void stopForceMotor();
void debugPrintJoystick(const int16_t diffX, const int16_t diffY);

void runLedClockwise();
void runLedClockwiseContinuous();
void runLedCounterclockwise();
void runLedCounterclockwiseContinuous();
void runLedBlink();

// -----------------------------------------------------------------------------
// Timer helpers
// -----------------------------------------------------------------------------
uint16_t calculateTimerTop(const uint32_t fosc,
                           const uint16_t prescaler,
                           const uint32_t timeTopUs) {
  const uint32_t scaled = ((fosc / 1000000UL) * timeTopUs) + (prescaler / 2UL);
  return static_cast<uint16_t>((scaled / prescaler) - 1UL);
}

void initTimer5(const uint32_t microSeconds) {
#if !defined(__AVR_ATmega2560__)
#error Only ATMEGA2560 MCU supported (Arduino Mega 2560)
#endif

  uint8_t prescaleBits = 0;
  uint32_t timerTicks = ((F_CPU / 1000000UL) * microSeconds) - 1UL;

  if (timerTicks < kTimerMax) {
    prescaleBits = (1 << CS50);
    timerTicks = calculateTimerTop(F_CPU, 1, microSeconds);
  } else if ((timerTicks >>= 3) < kTimerMax) {
    prescaleBits = (1 << CS51);
    timerTicks = calculateTimerTop(F_CPU, 8, microSeconds);
  } else if ((timerTicks >>= 3) < kTimerMax) {
    prescaleBits = (1 << CS51) | (1 << CS50);
    timerTicks = calculateTimerTop(F_CPU, 64, microSeconds);
  } else if ((timerTicks >>= 2) < kTimerMax) {
    prescaleBits = (1 << CS52);
    timerTicks = calculateTimerTop(F_CPU, 256, microSeconds);
  } else if ((timerTicks >>= 2) < kTimerMax) {
    prescaleBits = (1 << CS52) | (1 << CS50);
    timerTicks = calculateTimerTop(F_CPU, 1024, microSeconds);
  } else {
    prescaleBits = (1 << CS52) | (1 << CS50);
    timerTicks = kTimerMax;
  }

  TCCR5B &= ~((1 << CS52) | (1 << CS51) | (1 << CS50));
  TCCR5A = 0x00;
  TCCR5B = (1 << WGM52);
  TCCR5C = 0x00;
  OCR5A = static_cast<uint16_t>(timerTicks);
  TIMSK5 = (1 << OCIE5A);
  TCNT5 = 0x0000;
  TCCR5B |= prescaleBits;
}

ISR(TIMER5_COMPA_vect) {
  ADCSRA |= (1 << ADSC);
}

// -----------------------------------------------------------------------------
// ADC helpers
// -----------------------------------------------------------------------------
void initADC() {
#if !defined(__AVR_ATmega2560__)
#error Only ATMEGA2560 MCU supported (Arduino Mega 2560)
#endif

  uint8_t tempCh = (gCurrentAdcChannel - 54) & 0x0F;

  ADCSRA = 0x00;
  ADMUX = 0x00;
  ADMUX |= (1 << REFS0);
  ADMUX |= (tempCh & ((1 << MUX2) | (1 << MUX1) | (1 << MUX0)));

  ADCSRB = 0x00;
  ADCSRB |= (tempCh & (1 << MUX5));

  ADCSRA |= (1 << ADEN) | (1 << ADIE) | (1 << ADPS2) | (1 << ADPS1) | (1 << ADPS0);
}

ISR(ADC_vect) {
  const uint16_t tempAdc = ADC;
  uint8_t tempCh;

  ADMUX &= ~((1 << MUX4) | (1 << MUX3) | (1 << MUX2) | (1 << MUX1) | (1 << MUX0));
  ADCSRB &= ~(1 << MUX5);

  if (gCurrentAdcChannel == kPosXPin) {
    gPosX = static_cast<int16_t>(tempAdc);
    gCurrentAdcChannel = kPosYPin;
    tempCh = (gCurrentAdcChannel - 54) & 0x0F;
    ADMUX |= (tempCh & ((1 << MUX2) | (1 << MUX1) | (1 << MUX0)));
    ADCSRB |= (tempCh & (1 << MUX5));
    ADCSRA |= (1 << ADSC);
  } else if (gCurrentAdcChannel == kPosYPin) {
    gPosY = static_cast<int16_t>(tempAdc);
    gCurrentAdcChannel = kPosXPin;
    tempCh = (gCurrentAdcChannel - 54) & 0x0F;
    ADMUX |= (tempCh & ((1 << MUX2) | (1 << MUX1) | (1 << MUX0)));
    ADCSRB |= (tempCh & (1 << MUX5));
    gAdcSampleReady = true;
  }
}

// -----------------------------------------------------------------------------
// External interrupt
// -----------------------------------------------------------------------------
void handleSlideTrigger() {
  gSlideResetRequested = true;
}

// -----------------------------------------------------------------------------
// Arduino lifecycle
// -----------------------------------------------------------------------------
void setup() {
  if (kEnableSerialDebug) {
    Serial.begin(kSystemBaudRate);
  }
  Serial1.begin(kBpodBaudRate);

  gMotorShield.begin(MOTOR_AB);
  gMotorShield.set_forward_direction(MOTOR_A, LOW);
  gMotorShield.set_forward_direction(MOTOR_B, LOW);

  pinMode(kForceEnablePin, INPUT);
  pinMode(kForceOnsetMonitorPin, INPUT);
  pinMode(kSlideTrigPin, INPUT);

  pinMode(kLedDataOutPin, OUTPUT);
  pinMode(kLedSeq1TrigPin, INPUT);
  pinMode(kLedSeq2TrigPin, INPUT);
  pinMode(kInternalLedPin, OUTPUT);

  pinMode(kBuzzerTrigPin, INPUT);
  pinMode(kBuzzerOutPin, OUTPUT);

  pinMode(kTtlReachPin, OUTPUT);
  pinMode(kTtlForceOnPin, OUTPUT);
  pinMode(kTtlTooEarlyPin, OUTPUT);

  digitalWrite(kTtlReachPin, LOW);
  digitalWrite(kTtlForceOnPin, LOW);
  digitalWrite(kTtlTooEarlyPin, LOW);

  FastLED.addLeds<WS2811, kLedDataOutPin, RGB>(gLeds, kNumLeds);
  FastLED.setBrightness(kLedBrightness);
  FastLED.clear(true);

  attachInterrupt(digitalPinToInterrupt(kSlideTrigPin), handleSlideTrigger, RISING);

  gCurrentAdcChannel = kPosXPin;
  initADC();
  initTimer5(kAdcPeriodUs);
}

void loop() {
  if (gSlideResetRequested) {
    runSlideCloseSequence();
    gSlideResetRequested = false;
  }

  readBpodForceCommand();

  if (gAdcSampleReady) {
    processLatestJoystickSample();
  }

  updateLedSequences();
}

// -----------------------------------------------------------------------------
// Main task logic
// -----------------------------------------------------------------------------
void readBpodForceCommand() {
  if (Serial1.available() > 0) {
    gForceFactor = static_cast<uint8_t>(Serial1.read());
  }
}

void processLatestJoystickSample() {
  int16_t posX = 0;
  int16_t posY = 0;

  noInterrupts();
  gAdcSampleReady = false;
  posX = gPosX;
  posY = gPosY;
  interrupts();

  if (gBaselineResetRequested) {
    recalibrateBaseline();
    gBaselineResetRequested = false;
  }

  processJoystickSample(posX, posY);
}

void processJoystickSample(const int16_t posX, const int16_t posY) {
  noInterrupts();
  const int16_t baselineX = gBaselineX;
  const int16_t baselineY = gBaselineY;
  interrupts();

  gDiffX = abs(posX - baselineX);
  gDiffY = posY - baselineY;

  debugPrintJoystick(gDiffX, gDiffY);

  updateForceMotor(gDiffX, gDiffY);
  updateTooEarlyPulse(gDiffX, gDiffY);
  updateBuzzerAndReachTtl(gDiffX, gDiffY);
}

void updateForceMotor(const int16_t /*diffX*/, const int16_t diffY) {
  if (diffY > kForceOnsetThresholdY) {
    if (digitalRead(kForceOnsetMonitorPin) == HIGH) {
      // Return repeated brief pulses to Bpod while the onset-monitor gate is
      // high and the joystick remains above threshold. Repeated pulses are
      // intentional here to reduce the chance that Bpod misses a single pulse.
      sendPulse(kTtlForceOnPin, kTtlPulseMs);
    }

    if (digitalRead(kForceEnablePin) == HIGH) {
      gMotorShield.backward(MOTOR_A, gForceFactor);
      delayMicroseconds(kMotorCommandPauseUs);
      return;
    }
  }

  stopForceMotor();
}

void updateTooEarlyPulse(const int16_t diffX, const int16_t diffY) {
  if ((diffY > kTooEarlyThresholdY) && (diffX < kXWindowLimit)) {
    // Return repeated brief pulses to Bpod while the too-early
    // criterion remains true. Repeated pulses are intentional here.
    sendPulse(kTtlTooEarlyPin, kTtlPulseMs);
  }
}

void updateBuzzerAndReachTtl(const int16_t diffX, const int16_t diffY) {
  // kBuzzerTrigPin is a Bpod-controlled state gate. When HIGH, Arduino
  // enters the active monitoring state: it evaluates joystick position,
  // drives position-dependent buzzer output, and returns reach detection
  // to Bpod on kTtlReachPin.
  if (digitalRead(kBuzzerTrigPin) == HIGH) {
    if ((diffX < kBuzzerXGate) && (diffY > kBuzzerYGate)) {
      tone(kBuzzerOutPin, static_cast<unsigned int>((diffY * diffY) / 5 + kBaseBuzzerHz));
    } else {
      tone(kBuzzerOutPin, kBaseBuzzerHz);
    }

    // Return a sustained HIGH signal to Bpod while the reach criterion is
    // satisfied during the active monitoring state.
    if ((diffY > kForceOffsetThresholdY) && (diffX < kXWindowLimit)) {
      digitalWrite(kTtlReachPin, HIGH);
    } else {
      digitalWrite(kTtlReachPin, LOW);
    }
  } else {
    // Outside the Bpod monitoring state, buzzer output is disabled and the
    // reach-return line is held LOW.
    noTone(kBuzzerOutPin);
    digitalWrite(kTtlReachPin, LOW);
  }
}

void updateLedSequences() {
  if (digitalRead(kLedSeq1TrigPin) == HIGH) {
    runLedCounterclockwise();
  }

  if (digitalRead(kLedSeq2TrigPin) == HIGH) {
    runLedBlink();
  }
}

void recalibrateBaseline() {
  // Re-estimate joystick resting baseline during the Bpod dummy/reset state
  // using a fresh average of new samples only.
  delay(2000);

  int32_t sumX = 0;
  int32_t sumY = 0;

  for (uint8_t ii = 0; ii < 10; ++ii) {
    noInterrupts();
    const int16_t posX = gPosX;
    const int16_t posY = gPosY;
    interrupts();

    sumX += posX;
    sumY += posY;
    delay(200);
  }

  noInterrupts();
  gBaselineX = static_cast<int16_t>(sumX / 10);
  gBaselineY = static_cast<int16_t>(sumY / 10);
  interrupts();
}

void runSlideCloseSequence() {
  gMotorShield.forward(MOTOR_B, 255);
  delay(240);
  gMotorShield.stop(MOTOR_B);
  delay(500);

  gMotorShield.backward(MOTOR_B, 255);
  delay(230);
  gMotorShield.stop(MOTOR_B);
  delayMicroseconds(10);

  digitalWrite(kTtlReachPin, LOW);
  gBaselineResetRequested = true;
}

void sendPulse(const uint8_t pin, const uint16_t pulseMs) {
  digitalWrite(pin, HIGH);
  delay(pulseMs);
  digitalWrite(pin, LOW);
}

void stopForceMotor() {
  gMotorShield.stop(MOTOR_A);
  delayMicroseconds(kMotorCommandPauseUs);
}

void debugPrintJoystick(const int16_t diffX, const int16_t diffY) {
  if (!kEnableSerialDebug) {
    return;
  }

  static uint32_t lastPrintMs = 0;
  const uint32_t now = millis();
  if (now - lastPrintMs < 10) {
    return;
  }
  lastPrintMs = now;

  Serial.print(diffY);
  Serial.print(' ');
  Serial.println(diffX);
}

// -----------------------------------------------------------------------------
// LED helper routines (kept largely unchanged from v2)
// -----------------------------------------------------------------------------
void runLedClockwise() {
  for (int whiteLed = 0; whiteLed < kNumLeds - 4; ++whiteLed) {
    gLeds[whiteLed] = CRGB::White;
    gLeds[whiteLed + 1] = CRGB::White;
    gLeds[whiteLed + 2] = CRGB::White;
    gLeds[whiteLed + 3] = CRGB::White;
    gLeds[whiteLed + 4] = CRGB::White;
    FastLED.show();

    gLeds[whiteLed] = CRGB::Black;
    gLeds[whiteLed + 1] = CRGB::Black;
    gLeds[whiteLed + 2] = CRGB::Black;
    gLeds[whiteLed + 3] = CRGB::Black;
    gLeds[whiteLed + 4] = CRGB::Black;
    FastLED.show();
  }
}

void runLedClockwiseContinuous() {
  gLeds[0] = CRGB::White;
  gLeds[1] = CRGB::White;
  gLeds[2] = CRGB::White;
  gLeds[3] = CRGB::White;
  gLeds[4] = CRGB::White;
  FastLED.show();
  delay(0);

  for (int ledCntr = 1; ledCntr < kNumLeds - 4; ++ledCntr) {
    gLeds[ledCntr - 1] = CRGB::Black;
    gLeds[ledCntr + 4] = CRGB::White;
    FastLED.show();
    delay(0);
  }
}

void runLedCounterclockwise() {
  for (int whiteLed = kNumLeds - 5; whiteLed >= 0; --whiteLed) {
    gLeds[whiteLed] = CRGB::White;
    gLeds[whiteLed + 1] = CRGB::White;
    gLeds[whiteLed + 2] = CRGB::White;
    gLeds[whiteLed + 3] = CRGB::White;
    gLeds[whiteLed + 4] = CRGB::White;
    FastLED.show();
    delay(5);

    gLeds[whiteLed] = CRGB::Black;
    gLeds[whiteLed + 1] = CRGB::Black;
    gLeds[whiteLed + 2] = CRGB::Black;
    gLeds[whiteLed + 3] = CRGB::Black;
    gLeds[whiteLed + 4] = CRGB::Black;
    FastLED.show();
    delay(5);
  }
}

void runLedCounterclockwiseContinuous() {
  gLeds[(kNumLeds - 1) - 0] = CRGB::White;
  gLeds[(kNumLeds - 1) - 1] = CRGB::White;
  gLeds[(kNumLeds - 1) - 2] = CRGB::White;
  gLeds[(kNumLeds - 1) - 3] = CRGB::White;
  gLeds[(kNumLeds - 1) - 4] = CRGB::White;
  FastLED.show();
  delay(0);

  for (int ledCntr = kNumLeds - 2; ledCntr >= 4; --ledCntr) {
    gLeds[ledCntr + 1] = CRGB::Black;
    gLeds[ledCntr - 4] = CRGB::White;
    FastLED.show();
    delay(0);
  }
}

void runLedBlink() {
  FastLED.setBrightness(2);
  fill_solid(gLeds, kNumLeds, CRGB::White);
  FastLED.show();
  delay(2);

  fill_solid(gLeds, kNumLeds, CRGB::Black);
  FastLED.show();
  delay(50);

  FastLED.setBrightness(kLedBrightness);
}