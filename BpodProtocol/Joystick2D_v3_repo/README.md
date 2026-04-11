# Joystick2D_v3 Bpod protocol

MATLAB/Bpod protocol for 2D joystick reaching with optional force-field perturbations, online NI-DAQ photometry acquisition, and live behavior/trajectory/photometry plots.

## Task overview
In the baseline null-field phase, each trial follows:

**4 s pre-cue period -> 0.5 s cue -> 1 s delay -> 5.5 s response window**

Mice must move the joystick from the rest zone to the target zone and hold there for the configured hold time (`S.GUI.HoldingTime`; default `0.15 s` in this package). Successful trials enter a 1 s pre-reward delay followed by water delivery. Early reaches before the go period trigger a warning buzzer and terminate the trial. The between-trial interval consists of a fixed 2 s break plus a randomized truncated-exponential component between 1 and 6 s (total 3-8 s).

`Joystick2D_v3_Phase.m` defines additional phases for cued force, unexpected force, reward omission, optogenetic, and left/right force conditions.

## File structure
- `Joystick2D_v3.m` — main protocol and Bpod state machine
- `Bpod_Joystick2D_v3_parameter.m` — GUI defaults
- `Joystick2D_v3_Phase.m` — phase-specific trial tables
- `BpodParam_PCdep_Joystick_v3.m` — computer-specific NI-DAQ/LED defaults
- `Nidaq_*.m` — NI-DAQ acquisition and LED modulation helpers
- `Online_*.m` — online lick, photometry, and trajectory plots
- `GenerateBeatWave*.m`, `WhiteNoise.m`, `beepsound.m` — sound generation
- `SoftCodeHandler_PlaySound.m` — PsychToolbox sound playback

## Run order
1. Open Bpod and place this folder on the MATLAB path.
2. Edit `BpodParam_PCdep_Joystick_v3.m` so the computer name, NI device, and LED amplitudes match your rig.
3. Launch `Joystick2D_v3`.
4. Choose the task phase and GUI settings, then start the session.
5. The protocol initializes sounds, generates a weighted trial sequence, starts NI-DAQ background acquisition on each trial, runs the Bpod state matrix, saves trial data, and updates the online figures.

## Dependencies
- Bpod MATLAB framework
- PsychToolbox (`PsychToolboxSoundServer`)
- MATLAB Data Acquisition Toolbox with NI support
- National Instruments DAQ hardware
- Optional: PulsePal for optogenetic sessions (`S.GUI.OptoStim = 1`)

## Notes
- Cue/go hardware can be delivered externally; this MATLAB code sends the timing outputs and trial-state logic.
- If you use optogenetic sessions, update the hard-coded PulsePal COM port and `.mat` file path in `Joystick2D_v3.m`.
- Review the `WireState`, `Wire1High`, `Wire2High`, and `Wire3High` mappings against your joystick/rest-target hardware before running.
