# Joystick Force Field Controller (Arduino Mega 2560)

This repository contains the Arduino control code for a joystick-based force-field behavioral apparatus used with **Bpod**. The sketch runs on an **Arduino Mega 2560** and coordinates joystick sampling, motor control, return TTL signals to Bpod, position-dependent buzzer output, and LED feedback.

The current release file is:

- `Joystick_ForceField_v3_ready_to_run.ino`

## What the code does

The sketch performs five core functions:

1. Samples joystick position (`X`, `Y`) at **1 kHz** using **Timer5** and ADC interrupts.
2. Receives a one-byte motor amplitude command from **Bpod** over `Serial1`.
3. Drives **Motor A** continuously during Bpod-defined force states when joystick displacement exceeds the force-on threshold.
4. Returns TTL signals to **Bpod** for force-on detection, too-early movement detection, and successful reach detection.
5. Drives a **position-dependent buzzer** and **LED sequences** during Bpod-defined task states.

## Behavioral logic

### 1) Joystick sampling and baseline

- Joystick analog signals are read from `A8` (`X`) and `A9` (`Y`).
- The code stores a resting joystick baseline as `gBaselineX` and `gBaselineY`.
- Displacement is computed as:
  - `diffX = abs(posX - baselineX)`
  - `diffY = posY - baselineY`

A baseline reset is triggered after the slide-reset sequence. During this reset, the sketch waits for the joystick to settle and then computes a **fresh average of 10 new samples** to define the new resting baseline.

### 2) Force motor control

**Motor A** is the force-field motor.

Force is applied when both of the following are true:

- `diffY > kForceOnsetThresholdY`
- `kForceEnablePin` is `HIGH`

When those conditions are satisfied, the sketch drives **Motor A backward** using the force amplitude most recently received from Bpod over `Serial1`.

When either condition is not satisfied, Motor A is stopped.

### 3) Force-on return signal to Bpod

`kTtlForceOnPin` is a return line to Bpod.

When:

- `diffY > kForceOnsetThresholdY`, and
- `kForceOnsetMonitorPin` is `HIGH`

Arduino sends **repeated 2 ms TTL pulses** back to Bpod while the condition remains true.

This repeated pulsing is intentional. It is used as a robustness measure so the Bpod state machine is less likely to miss a single brief edge event.

### 4) Too-early return signal to Bpod

`kTtlTooEarlyPin` is another return line to Bpod.

When:

- `diffY > kTooEarlyThresholdY`, and
- `diffX < kXWindowLimit`

Arduino sends **repeated 2 ms TTL pulses** back to Bpod while the condition remains true.

As above, repeated pulsing is intentional and is used to reduce the chance of missed event detection by Bpod.

### 5) Active monitoring state, buzzer output, and reach detection

`kBuzzerTrigPin` is a **Bpod-controlled state gate**.

When `kBuzzerTrigPin` is `HIGH`, Arduino enters the **active monitoring state**. During this state, the sketch does two things at the same time:

- evaluates joystick position for successful reach detection
- drives buzzer output based on joystick position

#### Buzzer logic

During the active monitoring state:

- if `diffX < kBuzzerXGate` and `diffY > kBuzzerYGate`, the buzzer frequency is updated as a function of joystick position:
  - `frequency = (diffY * diffY) / 5 + kBaseBuzzerHz`
- otherwise, the buzzer plays a baseline tone at `kBaseBuzzerHz`

When `kBuzzerTrigPin` is `LOW`, the buzzer is turned off.

#### Reach return line

During the same active monitoring state, Arduino evaluates the reach criterion:

- `diffY > kForceOffsetThresholdY`
- `diffX < kXWindowLimit`

If that criterion is satisfied, Arduino sets `kTtlReachPin` **HIGH**.
If the criterion is not satisfied, Arduino sets `kTtlReachPin` **LOW**.

This is a **sustained level signal**, not a pulse.

The intended use is that Bpod enters a state in which it monitors actual joystick position, and Arduino returns a HIGH level on `kTtlReachPin` when the joystick enters the valid reach zone, allowing Bpod to transition to the next state.

When `kBuzzerTrigPin` is `LOW`, `kTtlReachPin` is forced `LOW`.

### 6) Slide reset and baseline recalibration

`kSlideTrigPin` is connected to a Bpod output and is configured as an interrupt input.

When a rising edge is detected:

- **Motor B** runs forward briefly
- Motor B stops
- Motor B runs backward briefly
- Motor B stops again
- `kTtlReachPin` is set LOW
- baseline recalibration is requested

This sequence is used to recenter/reset the joystick apparatus before the next trial block or monitoring period.

### 7) LED outputs

Two Bpod-controlled digital inputs trigger LED routines:

- `kLedSeq1TrigPin` triggers a counterclockwise LED sequence
- `kLedSeq2TrigPin` triggers a blink sequence

These LED helper routines are carried over largely unchanged from the working task code.

## Hardware and dependencies

### Board

- **Arduino Mega 2560**

### Libraries

The sketch depends on:

- `FastLED`
- `src/MotorShield/MotorShield.h`

If your setup also uses additional local support libraries, keep them in the repository or provide installation instructions.

## Pin map

### Inputs from Bpod / hardware

- `kLedSeq1TrigPin = 22` — Bpod output for LED sequence 1
- `kLedSeq2TrigPin = 23` — Bpod output for LED sequence 2
- `kSlideTrigPin = 2` — Bpod output for slide reset trigger
- `kBuzzerTrigPin = 24` — Bpod state gate for active monitoring and buzzer output
- `kForceOnsetMonitorPin = 25` — Bpod gate for force-on return pulses
- `kForceEnablePin = 26` — Bpod gate enabling continuous Motor A force output
- `kPosXPin = A8` — joystick X analog input
- `kPosYPin = A9` — joystick Y analog input

### Outputs from Arduino

- `kBuzzerOutPin = 36` — buzzer output
- `kTtlReachPin = 30` — sustained reach return line to Bpod
- `kTtlForceOnPin = 31` — repeated force-on pulse return line to Bpod
- `kTtlTooEarlyPin = 32` — repeated too-early pulse return line to Bpod
- `kLedDataOutPin = 37` — WS2811 LED data output
- `kInternalLedPin = 13` — onboard status LED pin

### Motors

- **Motor A** — force-field motor
- **Motor B** — slide-close / reset motor

## Key parameters

Default values in the current sketch:

- `kAdcPeriodUs = 1000` → 1 kHz sampling
- `kTtlPulseMs = 2` → Bpod return-pulse width
- `kForceOnsetThresholdY = 20`
- `kForceOffsetThresholdY = 50`
- `kTooEarlyThresholdY = 100`
- `kXWindowLimit = 60`
- `kBuzzerXGate = 50`
- `kBuzzerYGate = 40`
- `kBaseBuzzerHz = 2000`
- `kMotorCommandPauseUs = 300`

These values are in ADC units unless otherwise noted.

## Serial communication

The sketch uses two serial interfaces:

- `Serial1` at **115200 baud** for incoming Bpod force commands
- `Serial` at **115200 baud** for optional USB debug output / Serial Plotter

In the current release, debug output is enabled with:

```cpp
constexpr bool kEnableSerialDebug = true;
```

The sketch prints two values to USB serial for plotting:

- first column: `diffY`
- second column: `diffX`

The debug stream is throttled to approximately **100 Hz** to improve compatibility with Arduino Serial Plotter.

## Setup and upload

1. Install the required libraries.
2. Place the sketch in a folder named exactly like the `.ino` file.
3. Open the sketch in the Arduino IDE.
4. Select:
   - **Board:** Arduino/Genuino Mega or Mega 2560
   - **Processor:** ATmega2560 (Mega 2560)
5. Upload the sketch.
6. If needed, open **Tools → Serial Plotter** at **115200 baud** to inspect `diffY` and `diffX`.

## Suggested repository contents

A minimal repository for publication should include:

- `Joystick_ForceField_v3_ready_to_run.ino`
- `src/MotorShield/`
- `README.md`
- `LICENSE`
- a brief wiring diagram or `hardware_wiring.md`

## Notes for reuse

This code is written for a specific joystick + motor + Bpod behavioral rig. Anyone reusing it will likely need to adapt:

- joystick baseline values
- threshold values
- Bpod state timing
- motor directions
- motor power scaling
- wiring and pin assignments

## Code availability statement template

Custom Arduino code used for joystick position sampling, motorized force-field control, TTL communication with Bpod, position-dependent buzzer output, and LED control is available in this repository. The repository includes the exact Arduino sketch used for behavioral control, along with required dependencies and instructions for compilation and upload.
