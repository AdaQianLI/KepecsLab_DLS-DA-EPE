# Joystick Force Field Controller

Arduino control code for a joystick-based force-field behavioral rig integrated with **Bpod**.

The sketch runs on an **Arduino Mega 2560** and performs four main jobs:

1. samples joystick position at **1 kHz**,
2. drives the force motor according to joystick displacement and Bpod state signals,
3. returns behavioral signals back to Bpod,
4. controls a go_cue output and LED routines during defined task states.

## Main modules

### 1. Joystick sampling and baseline
- Reads joystick analog inputs from `A8` (X) and `A9` (Y).
- Computes displacement relative to a stored resting baseline.
- Recalculates baseline after the slide-reset sequence using a fresh average of new samples.

### 2. Force motor control
- **Motor A** is the force-field motor.
- Bpod sends a one-byte motor amplitude command over `Serial1`.
- When joystick displacement exceeds the force-on threshold and the Bpod force-enable line is high, Arduino drives Motor A continuously.
- Otherwise, Motor A is stopped.

### 3. Return signals to Bpod
Arduino sends three behavioral signals back to Bpod:

- **Force-on return**: repeated brief TTL pulses while the force-on criterion is satisfied and the force-on monitor gate is active.
- **Too-early return**: repeated brief TTL pulses while the too-early movement criterion (`diffY > kTooEarlyThresholdY`) is satisfied.
- **Reach return**: a sustained HIGH level while the joystick is within the valid reach zone during the active monitoring state.

The repeated pulsing for force-on and too-early detection is intentional and helps reduce missed detection by the Bpod state machine.

### 4. Active monitoring state and go_cue
- `kBuzzerTrigPin` is a **Bpod-controlled state gate**.
- When this pin is high, Arduino enters the active monitoring state.
- During this state, Arduino:
  - evaluates joystick position for reach detection,
  - outputs a **go_cue** on `kBuzzerOutPin`.

The go_cue is state-dependent. Its **position-based frequency modulation is optional**:
- within a defined position window, cue frequency can be modulated as a function of joystick position,
- otherwise the cue remains at a baseline tone.

When the active monitoring gate is low, the go_cue is off and the reach return line is held low.

### 5. Slide reset
- A rising edge on `kSlideTrigPin` triggers the slide-reset sequence using **Motor B**.
- After reset, the sketch requests baseline recalibration before the next monitoring period.

### 6. LED routines
- Two Bpod-controlled input lines trigger LED routines.
- These routines provide task-related visual feedback and are kept as helper modules in the current release.

## Hardware
- **Board:** Arduino Mega 2560
- **Libraries:** `FastLED`, `src/MotorShield/MotorShield.h`

## Key I/O lines

### Inputs from Bpod
- `kBuzzerTrigPin` — active monitoring state gate and go_cue enable
- `kForceOnsetMonitorPin` — gate for force-on return pulses
- `kForceEnablePin` — enables continuous force output
- `kSlideTrigPin` — triggers slide reset
- `kLedSeq1TrigPin`, `kLedSeq2TrigPin` — LED routine triggers

### Outputs from Arduino
- `kTtlForceOnPin` — repeated force-on return pulses
- `kTtlTooEarlyPin` — repeated too-early return pulses
- `kTtlReachPin` — sustained reach return line
- `kBuzzerOutPin` — go_cue output
- `kLedDataOutPin` — LED data output

## Serial communication
- `Serial1` at **115200 baud** receives the motor amplitude command from Bpod.
- `Serial` at **115200 baud** can be used for USB debugging / Serial Plotter.

## Repository contents
A minimal release should include:
- `Joystick_ForceField_v3_ready_to_run.ino`
- `src/MotorShield/`
- `README.md`
- `LICENSE`
- optional wiring notes or diagram

## Code availability statement
Custom Arduino code used for joystick sampling, force-field motor control, TTL communication with Bpod, go_cue output, optional position-dependent cue modulation, and LED control is available in this repository.
