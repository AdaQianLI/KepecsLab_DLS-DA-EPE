#include <Arduino.h>
#include <FastLED.h>
#include "src/MotorShield/MotorShield.h"

/*******************************************************************************
 ***                             Define Constants                            ***
 ******************************************************************************/
#define SYSTEM_BAUDRATE 9600UL
#define BPOD_BAUDRATE 115200UL
#define MOTOR_A_RUN_DURATION_MS 4UL
#define MOTOR_A_PAUSE_DURATION_US 300UL
#define MOTOR_A_PAUSE_DURATION_US 300UL
#define X0 521
#define Y0 521
#define xlim 50
#define forceOffset 110  // forceOFF POS, wire1
#define forceOnset 40 // forceON POS, wire2
#define tooEalyPos 80 // 1st tooEarlycheck POS, wire3
#define NUM_LEDS 72
#define LED_BRIGHTNESS 10

/*******************************************************************************
 ***                             Define GPIOs                                ***
 ******************************************************************************/
#define LED_DATA_OUT_PIN 37  // Data pin for LEDs
#define LED_SEQ1_TRIG_PIN 22 // Bpod BNC OUT 1
#define LED_SEQ2_TRIG_PIN 23 // Bpod BNC OUT 2
#define SLIDE_TRIG_PIN 2     // Bpod Wire Out 1
#define BUZZER_TRIG_PIN 24   // Bpod Wire Out 2
#define interOnsetPin 25     // Bpod Wire Out 3
#define FORCE_ENABLE_PIN 26  // Bpod Wire Out 4
#define BUZZER_OUT_PIN 36    // Buzzer output
#define TTL 30               // Bpod Wire In 1
#define TTL_forceON 31       // Bpod Wire In 2
#define TTL_expect 32        // Bpod Wire In 3
#define POS_X_PIN A8         // Analog In 0 for x-axis of joystick
#define POS_Y_PIN A9         // Analog In 1 for y-axis of joystick
#define INTERNAL_LED 13      // Pin of the Arduino's internal LED

/*******************************************************************************
 ***                   Define constants and flags etc                        ***
 ******************************************************************************/
#define TIMERMAX 0xFFFF            // Max value for 16-bit timer
#define TTL_FLAG_INTER_ONSET 0x01  // bitmap for TTL_forceON flag
#define TTL_FLAG_TTL_EXPECT  0x02  // bitmap for TTL_expect flag

/*******************************************************************************
 ***                         Function Prototypes                             ***
 ******************************************************************************/
void Slideopen();
void Slideclose();
void LEDclockwise();
void LEDclockwise2();
void LEDCounterclockwise();
void LEDCounterclockwise2();
void LEDblink();

/*******************************************************************************
 ***                       Define Global Variables                           ***
 ******************************************************************************/
// Speed (force) for the force-feedback
uint8_t forceFactor = 0;
uint8_t forceFactor2 = 60;
// x- and y-position of joystick
// - must be volatile now since we access them in the ISR and man loop!!!
volatile int16_t posX;
volatile int16_t posY;

// difference between current and zero position
int16_t diffX;
int16_t diffY;

// LED array
CRGB leds[NUM_LEDS];

// Trigger variable for closing the slide
// - needs to be volatile since being used in ISR and main loop
volatile uint8_t closeTrigger = 0;

// Variable storing the next(!!!) channel to be sampled!
// - needs to be volatile since being used in ISR and main loop
volatile uint8_t currentADCch;

// Trigger variable to indicate that new adc values are sampled
// - needs to be volatile since being used in ISR and main loop
volatile uint8_t adcTrigger = 0;

// Timestamps for creating TTL signals
uint32_t interOnsetTimestamp;
uint32_t ttlExpectTimestamp;

// Flag variable for creating TTL pulses
uint8_t TTLflags;

// Create an object to access the motor driver board
MotorShield motorShield;

// Calculate the value for the compare registers or to preload tcnt...
uint16_t calculateTimerTop(const uint32_t fosc, 
                           const uint16_t prescaler,
                           const uint32_t timeTop)
{
  uint16_t returnValue;
  // Calculation with "intelligent" rounding
  returnValue = ((((fosc/1000000) * timeTop) + (prescaler/2))/(prescaler))-1;
  return (returnValue);
}

// Initialization for Timer5 for CTC mode (Clear timer on Compare Match)
void initTimer5(const uint32_t microSeconds)
{
  #if !defined(__AVR_ATmega2560__)
  #error Only ATMEGA2560 MCU supported (Ardiuno MEGA)
  #endif

  uint8_t prescale;
  uint32_t timerTicks;

  // Calculate how many timer ticks (without prescaler!) would be necessary to 
  // count to desired number of microseconds
  timerTicks = ((F_CPU/1000000) * microSeconds)-1;

  if (timerTicks < TIMERMAX)
  {
    // the value fits into 16-bits ->  prescaler to 1
    prescale = (1 << CS50);
    // calculate the actual value for OCR register
    timerTicks = calculateTimerTop(F_CPU, 1, microSeconds);
  }
  else if ((timerTicks >>= 3) < TIMERMAX) // Shift 3 pos right -> divide by 8
  {
    // the value fits into 16-bits if prescaler set to 8
    prescale = (1 << CS51);
    // calculate the actual value for OCR register
    timerTicks = calculateTimerTop(F_CPU, 8, microSeconds);
  }
  else if ((timerTicks >>= 3) < TIMERMAX) // Shift another 3 pos -> divide by 64
  {
    // the value fits into 16-bits if prescaler set to 64
    prescale = ((1 << CS51) | (1 << CS50));
    // calculate the actual value for OCR register
    timerTicks = calculateTimerTop(F_CPU, 64, microSeconds);
  }
  else if ((timerTicks >>= 2) < TIMERMAX) // Shift another 2 pos -> divide by 256
  {
    // the value fits into 16-bits if prescaler set to 256
    prescale = (1 << CS52);
    // calculate the actual value for OCR register
    timerTicks = calculateTimerTop(F_CPU, 256, microSeconds);
  }
  else if ((timerTicks >>= 2) < TIMERMAX) // Shift another 2 pos -> divide by 1024
  {
    // the value fits into 16-bits if prescaler set to 1024
    prescale = ((1 << CS52) | (1 << CS50));
    // calculate the actual value for OCR register
    timerTicks = calculateTimerTop(F_CPU, 1024, microSeconds);
  }
  else
  {
    // even if using a prescaler of 1024 the value does not fit into 16-bits...
    // use max value (65536 ticks) with prescaler set to 1024...
    prescale = ((1 << CS52) | (1 << CS50));
    timerTicks = TIMERMAX;
  }
  

  // We are sing Timer5 in Clear Timer on Compare MAtch mode - CTC
  // This is mode 4 and the WGM bits in the control registers need to be:
  // WGM5[3:0] = 0b0100

  // First, stop Timer5 in case it was running before (should not have happened)
  // by setting the Clock Select bits in TCCR5B to 0
  TCCR5B &= ~((1 << CS52) | (1 << CS51) | (1 << CS50));

  // Set TCCR5A (Timer/Counter Control Register Timer5 A)
  // +--------+--------+--------+--------+--------+--------+-------+-------+
  // | COM5A1 | COM5A0 | COM5B1 | COM5B0 | COM5C1 | COM5C0 | WGM51 | WGM50 |
  // +--------+--------+--------+--------+--------+--------+-------+-------+
  //
  // Since we are not using the OC5A/OC5B/OC5C output pins, the first 6 bits
  // need to be set to 0
  // - COM5A[1:0]: 0b00 
  // - COM5B[1:0]: 0b00
  // - COM5C[1:0]: 0b00
  //
  // Furthermore, we are using Timer5 in CTC mode WGM5[3:0]: 0b0100
  // - WGM5[1:0]: 0b00
  //
  // So TCCR5A must be set to 0x00 = (0b00000000)
  TCCR5A = 0x00;

  // Set TCCR5B (Timer/Counter Control Register Timer5 B)
  // +-------+-------+---+-------+-------+------+------+------+
  // | ICNC5 | ICES5 | - | WGM53 | WGM52 | CS52 | CS51 | CS50 |
  // +-------+-------+---+-------+-------+------+------+------+
  // 
  // Since we are not using the input capture capabilities of Timer5, the first 
  // 2 bits need to be set to 0
  // - ICNC5: 0b0   (Input Capture Noise Canceler)
  // - ICES5: 0b0   (Input Capture Edge Select)
  //
  // Set the 2 most significant bit of the WGM5[3:0] set accordgingly
  // - WGM5[3:2]: 0b01 (mode 4 of the waveform generator mode bit structure)
  //
  // The bits CS5[2:0] specify the prescaler. When CS5[2:0] set to any value 
  // except 000 the timer starts - so set them at the end of this function!
  TCCR5B = (1 << WGM52);

  // Set TCCR5C (Timer/Counter Control Register Timer5 C)
  // +-------+-------+-------+---+---+---+---+---+
  // | FOC5A | FOC5B | FOC3C | - | - | - | - | - |
  // +-------+-------+-------+---+---+---+---+---+
  // 
  // Since we are not using the force output compare options here, we can set 
  // set the first 3 bits to 0
  // - FOC5A: 0b0   (Force Output Compare for Channel A)
  // - FOC5B: 0b0   (Force Output Compare for Channel B)
  // - FOC3C: 0b0   (Force Output Compare for Channel C
  //
  // Since the other bits in this register are unused, the resulting value for 
  // this register are 0b00000000 (0x00)
  TCCR5C = 0x00;

  // We are using Timer5 in CTC mode. That means whenever the timer reaches the 
  // value stored in OCR5A, the timer will be reset and starts again.
  OCR5A = timerTicks;

  // Enable the interrupt for the Output Compare Match A event for Timer5
  TIMSK5 = (1 << OCIE5A);

  // Reset Timer/Counter 5
  TCNT5 = 0x0000;

  // Set prescaler value and start the timer
  TCCR5B |= prescale;
}

// Interrupt Service Routine for Timer5 Compare Match A interrupt
ISR(TIMER5_COMPA_vect)
{
  // Only start the ADC here
  ADCSRA |= (1 << ADSC);
}

// Initialization for ADC - it is actually a re-initialization since the ADC was
// already initialized in the core Arduino functions
void initADC()
{
#if !defined(__AVR_ATmega2560__)
#error Only ATMEGA2560 MCU supported (Ardiuno MEGA)
#endif

  uint8_t tempCh = 0;

  // First of all, turn the ADC off by resetting ADCSRA completely
  ADCSRA = 0x00;

  // The base value for ADC0 specified in pins_Arduino.h is 54 and all other ADC
  // are following in an ascending order (ADC1 -> 55, ADC2 -> 56, ... ADC15 -> 69)
  tempCh = currentADCch - 54;

  // Make sure only the 4 least significant bits are being evaluated
  tempCh &= 0x0F;

  // Initialize the ADMUX register
  ADMUX = 0x00;          // Reset ADMUX register
                         // - Especially the ADLAR bit
  ADMUX |= (1 << REFS0); // Enable internal reference with cap at AREF

  // Set bits in ADMUX accordingly to the selected channel
  ADMUX |= (tempCh & ((1 << MUX2) | (1 << MUX1) | (1 << MUX0)));

  // Initialize the ADCSRB register
  ADCSRB = 0x00;                      // Reset the register
  ADCSRB |= (tempCh & ((1 << MUX5))); //

  // Finally configure the ADCSRA register
  ADCSRA |= ((1 << ADEN) | (1 << ADIE) | (1 << ADPS2) | (1 << ADPS1) | (1 << ADPS0));
  // ADEN:    Enable ADC (but don't start it yet)
  // ADIE:    Enable Interrupt
  // ADSP2:0: Set clock prescaler to 128 - ADC_clock = 125 kHz
}

// Interrupt Service Routine for ADC conversions
ISR(ADC_vect)
{
  uint16_t tempADC = 0;
  uint8_t tempCh;

  // Copy the ADC value into memory
  // Usually we'd have to first load ADCL and then ADCH but we can also directly
  // access both registers at the same time with the ADC SFR
  // (special function register) definition...
  tempADC = ADC;

  // Reset MUX bits in ADMUX and ADCSRB
  ADMUX  &= ~((1 << MUX4) | (1 << MUX3) | (1 << MUX2) | (1 << MUX1) | (1 << MUX0));
  ADCSRB &= ~(1 << MUX5);

  if (currentADCch == POS_X_PIN)
  {
    // Save ADC value for x position
    posX = tempADC;

    // Set currentADCch to pin associated with Y position of joystick
    currentADCch = POS_Y_PIN;

    // The base value for ADC0 specified in pins_Arduino.h is 54 and all other 
    // ADC channels are following in an ascending order
    // (ADC1 -> 55, ADC2 -> 56, ... ADC15 -> 69)
    tempCh = currentADCch - 54;

    // Make sure only the 4 least significant bits are being evaluated
    tempCh &= 0x0F;

    // Set bits in ADMUX and ADCSRB accordingly
    ADMUX  |= (tempCh & ((1 << MUX2) | (1 << MUX1) | (1 << MUX0)));
    ADCSRB |= (tempCh & ((1 << MUX5)));

    // Restart ADC - to get the Y position for this sampling interval
    ADCSRA |= (1 << ADSC);
  }

  else if (currentADCch == POS_Y_PIN)
  {
    // Save ADC value for y position
    posY = tempADC;

    // Set currentADCch to pin associated with Y position of joystick
    currentADCch = POS_X_PIN;

    // The base value for ADC0 specified in pins_Arduino.h is 54 and all other 
    // ADC channels are following in an ascending order
    // (ADC1 -> 55, ADC2 -> 56, ... ADC15 -> 69)
    tempCh = currentADCch - 54;

    // Make sure only the 4 least significant bits are being evaluated
    tempCh &= 0x0F;

    // Set bits in ADMUX and ADCSRB accordingly
    ADMUX |= (tempCh & ((1 << MUX2) | (1 << MUX1) | (1 << MUX0)));
    ADCSRB |= (tempCh & ((1 << MUX5)));

    // Do not restart ADC!
    // This will be done at Timer5's ISR to start a new cycle...

    // Set ADC trigger variable
    adcTrigger = 1;
  }
}

// Interrupt handler for the ext. interrupt event - Bpod triggering for slide
void interruptHandler()
{
  closeTrigger = 1;
}

// Initialization
void setup()
{
  // Initialize Serial interfaces
  Serial.begin(SYSTEM_BAUDRATE); // Serial interface for plotting on PC
  Serial1.begin(BPOD_BAUDRATE);  // Serial interface for communication with Bpod

  motorShield.begin(MOTOR_AB);
  motorShield.set_forward_direction(MOTOR_A, LOW);
  motorShield.set_forward_direction(MOTOR_B, LOW);

  // GPIOs for controlling force motor
  pinMode(FORCE_ENABLE_PIN, INPUT);

  // GPIOs for LEDs
  pinMode(LED_DATA_OUT_PIN, OUTPUT);
  pinMode(LED_SEQ1_TRIG_PIN, INPUT);
  pinMode(LED_SEQ2_TRIG_PIN, INPUT);
  pinMode(INTERNAL_LED, OUTPUT);

  // GPIOs for buzzer
  pinMode(BUZZER_TRIG_PIN, INPUT);
  pinMode(BUZZER_OUT_PIN, OUTPUT);

  pinMode(TTL, OUTPUT);
  pinMode(TTL_forceON, OUTPUT);
  pinMode(TTL_expect, OUTPUT);

  // Initialze LEDs
  FastLED.addLeds<WS2811, LED_DATA_OUT_PIN, RGB>(leds, NUM_LEDS);
  FastLED.setBrightness(LED_BRIGHTNESS);

  // Enable interrupt for centering joystick
  attachInterrupt(digitalPinToInterrupt(SLIDE_TRIG_PIN), interruptHandler, RISING);

  // Additions to measure time the joystick is above the threshold
  // ---------------------------------------------------------------------------
  // Set the current ADC channel to the x position input
  currentADCch = POS_X_PIN;
  // Re-initialize the ADC (it was already initialized in the core Arduino files)
  // Whenever the ADC finishes a conversion it raises an interrupt
  initADC();
  // Initialize and start Timer5 to generate an interrupt every 1 ms (1000 us). 
  // In the ISR a new ADC conversion will be started
  initTimer5(1000);
}

// The main "program"
void loop()
{
  // Check for external interupt event to close the slider
  if (closeTrigger == 1)
  {
    Slideclose();
    closeTrigger = 0;
    digitalWrite(TTL, LOW);
  }

  // ---------------------------------------------------------------------------
  // Do some time sensitive stuff here - START
  // ---------------------------------------------------------------------------

  // Reset TTL_forceON TTL pulse after 2 ms
 // if (TTLflags & TTL_FLAG_INTER_ONSET){
 //   if ( (micros() - interOnsetTimestamp ) > 5000)
 //   {
 //       digitalWrite(TTL_forceON, LOW);
 //       TTLflags &= ~TTL_FLAG_INTER_ONSET;
  //  }
  //}

  // Reset TTL_expect TTL pulse after 2 ms
//  if (TTLflags & TTL_FLAG_TTL_EXPECT){
 //   if ( (micros() - ttlExpectTimestamp ) > 5000)
 //   {
 //       digitalWrite(TTL_expect, LOW);
 //       TTLflags &= ~TTL_FLAG_TTL_EXPECT;
 //   }
 // }
  // ---------------------------------------------------------------------------
  // Do some time sensitive stuff here - END
  // ---------------------------------------------------------------------------

  // Check for incoming bytes from Bpod
  if (Serial1.available())
  {
    forceFactor = Serial1.read();
  }

  // All the commands in the following block will be triggered through the
  // 1ms timer
  if (adcTrigger)
  {
    // Reset ADC trigger variable
    adcTrigger = 0;

    // Convert joystick position in relation zero position
    diffX = abs(posX - X0);
    diffY = Y0 - posY;

    // Output values via serial terminal
   Serial.print(diffX);
   Serial.print(" ");
   Serial.println(diffY);


 
    if (diffY > forceOnset)
    {
      if (digitalRead(interOnsetPin) == HIGH)
      {
        // Send TTL pulse to Bpod
        digitalWrite(TTL_forceON, HIGH);
    //   interOnsetTimestamp = millis();
     //   TTLflags |= TTL_FLAG_INTER_ONSET;
        delay(2);
        digitalWrite(TTL_forceON, LOW);
      }
   // }
      if (digitalRead(FORCE_ENABLE_PIN) == HIGH)
      {
        motorShield.backward(MOTOR_A, forceFactor);
        // This delay can be deleted, right?!?
        delayMicroseconds(300);
      }
   // }
    else
    {
      motorShield.stop(MOTOR_A);
      // This delay can be deleted, right?!?
      delayMicroseconds(300);
    }
  }
      else
    {
      motorShield.stop(MOTOR_A);
      // This delay can be deleted, right?!?
      delayMicroseconds(300);
    }
          if ((diffY > tooEalyPos) && (diffX < xlim))
      {
        digitalWrite(TTL_expect, HIGH);
  //      ttlExpectTimestamp = millis();
   //     TTLflags |= TTL_FLAG_TTL_EXPECT;
        delay(2);
        digitalWrite(TTL_expect, LOW);
      }
    if (digitalRead(BUZZER_TRIG_PIN) == HIGH)
    {
      //
      if ((diffX < 50) && (diffY > 40))
      {
        tone(BUZZER_OUT_PIN, diffY * diffY / 5 + 2000); //*30
        // The following delay can be deleted since the outer code block will be 
        // evaluated every 1 ms anyhow... so we don't need to wait here...
        //delay(1);
      }
      else
      {
        tone(BUZZER_OUT_PIN, 2000);
        // The following delay can be deleted since the outer code block will be 
        // evaluated every 1 ms anyhow... so we don't need to wait here...
        //delay(1);
      }

      if ((diffY > forceOffset) && (diffX < xlim))
      {
        digitalWrite(TTL, HIGH);
      }
      else
      {
        digitalWrite(TTL, LOW);
      }
      

    }
    else
    {
      noTone(BUZZER_OUT_PIN);
      // The following delay can be deleted since the outer code block will be 
      // evaluated every 1 ms anyhow... so we don't need to wait here...
      //delay(1);
    }
  }

  // The following code is not affected by the ADC/Timer stuff so it needs to
  // be outside the if statement for new ADC values...

  if (digitalRead(LED_SEQ1_TRIG_PIN) == HIGH)
  {
    LEDCounterclockwise();
    //digitalWrite(INTERNAL_LED, HIGH);
  }
  else
  {
    //digitalWrite(INTERNAL_LED, LOW);
  }
  if (digitalRead(LED_SEQ2_TRIG_PIN) == HIGH)
  {
    LEDblink();
  }
}

void Slideclose()
{

  for (int i = 0; i <= 1200; i++)
  {
    motorShield.forward(MOTOR_B, 255);
    delayMicroseconds(160);
    motorShield.stop(MOTOR_B);
    delayMicroseconds(2);
  }
  delay(1000);
  //  motorShield.backward(MOTOR_B, 255);
  //  delay(340); //3400 * 0.1
  //  motorShield.stop(MOTOR_B);

  // for (int i = 0; i <= 6000; i++)
  // {
  //   motorShield.forward(MOTOR_B, 150);
  //   //delay(1);
  //   delayMicroseconds(80);
  // }
  for (int i = 0; i <= 1200; i++)
  {
    motorShield.backward(MOTOR_B, 255);
    delayMicroseconds(160); // 6000 * 0.08
    motorShield.stop(MOTOR_B);
    delayMicroseconds(2);
  }
}

void LEDclockwise()
{
  // Move a single white led
  for (int whiteLed = 0; whiteLed < NUM_LEDS - 4; whiteLed++)
  {
    // Turn our current led on to white, then show the leds
    leds[whiteLed] = CRGB::White;
    leds[whiteLed + 1] = CRGB::White;
    leds[whiteLed + 2] = CRGB::White;
    leds[whiteLed + 3] = CRGB::White;
    leds[whiteLed + 4] = CRGB::White;
    FastLED.show();
    //     delay(5);

    // Turn our current led back to black for the next loop around
    leds[whiteLed] = CRGB::Black;
    leds[whiteLed + 1] = CRGB::Black;
    leds[whiteLed + 2] = CRGB::Black;
    leds[whiteLed + 3] = CRGB::Black;
    leds[whiteLed + 4] = CRGB::Black;
    FastLED.show();
    //   delay(5);
  }
}

void LEDclockwise2()
{
  // Turn the first 5 LEDs on
  leds[0] = CRGB::White;
  leds[1] = CRGB::White;
  leds[2] = CRGB::White;
  leds[3] = CRGB::White;
  leds[4] = CRGB::White;
  FastLED.show();
  delay(0);

  // Move a single white led
  for (int ledCntr = 1; ledCntr < NUM_LEDS - 4; ledCntr++)
  {
    // Turn the first LED of the previous group off and the next LED on
    leds[ledCntr - 1] = CRGB::Black;
    leds[ledCntr + 4] = CRGB::White;
    FastLED.show();
    delay(0);
  }
}

void LEDCounterclockwise()
{
  for (int whiteLed = NUM_LEDS - 5; whiteLed >= 0; whiteLed--)
  {
    // Turn our current led on to white, then show the leds
    leds[whiteLed] = CRGB::White;
    leds[whiteLed + 1] = CRGB::White;
    leds[whiteLed + 2] = CRGB::White;
    leds[whiteLed + 3] = CRGB::White;
    leds[whiteLed + 4] = CRGB::White;
    FastLED.show();
    delay(5);

    // Turn our current led back to black for the next loop around
    leds[whiteLed] = CRGB::Black;
    leds[whiteLed + 1] = CRGB::Black;
    leds[whiteLed + 2] = CRGB::Black;
    leds[whiteLed + 3] = CRGB::Black;
    leds[whiteLed + 4] = CRGB::Black;
    FastLED.show();
    delay(5);
  }
}

void LEDCounterclockwise2()
{
  // Turn the last 5 leds onleds[whiteLed] = CRGB::White;
  leds[(NUM_LEDS - 1) - 0] = CRGB::White; // last LED
  leds[(NUM_LEDS - 1) - 1] = CRGB::White;
  leds[(NUM_LEDS - 1) - 2] = CRGB::White;
  leds[(NUM_LEDS - 1) - 3] = CRGB::White;
  leds[(NUM_LEDS - 1) - 4] = CRGB::White;
  FastLED.show();
  delay(0);

  for (int ledCntr = NUM_LEDS - 2; ledCntr >= 4; ledCntr--)
  {
    // Turn our current led on to white, then show the leds
    leds[ledCntr + 1] = CRGB::Black;
    leds[ledCntr - 4] = CRGB::White;
    FastLED.show();
    delay(0);
  }
}

void LEDblink()
{
  FastLED.setBrightness(2);
  // Turn all LEDs on simultaneously
  fill_solid(leds, NUM_LEDS, CRGB::White);
  FastLED.show();
  delay(2);

  // Turn all LEDs off simultaneously
  fill_solid(leds, NUM_LEDS, CRGB::Black);
  FastLED.show();
  delay(50);

  FastLED.setBrightness(LED_BRIGHTNESS);
}
