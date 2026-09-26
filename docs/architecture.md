# TCLab Project Architecture

## Control path

```text
Measured T1, T2
      │
      ▼
Output error ───────► Luenberger observer ─────► estimated state x̂
      │                                              │
      │                                              ▼
      └──────────── setpoint ───────► state-feedback controller K
                                                     │
                                                     ▼
                                              heater commands Q1,Q2
                                                     │
                                                     ▼
                                                   TCLab
```

The state-feedback controller uses the estimated state because the full internal state is not directly measured. The discrete Luenberger observer reconstructs the state from the measured temperatures and applied inputs.

## Identification path

```text
Measured input/output data
          │
          ├──────────────► MOESP ─────► linear state-space model
          │
          └──────────────► kernel regression ─────► nonlinear model
                                      │
                                      ▼
                              validation data
                                      │
                                      ▼
                         model / real-plant comparison
```

The identification workflow uses a linear MOESP model and a degree-2 polynomial-kernel model. Both are evaluated against held-out validation data, with an additional comparison against newly collected TCLab measurements from the original experiment.
