# Firefly task: perturbation, gain and density analyses

MATLAB analysis code behind

> Alefantis P., Lakshminarasimhan K., Avila E., Noel J.-P., Pitkow X., Angelaki D.E. (2022).
> **Sensory evidence accumulation using optic flow in a naturalistic navigation task.**
> *The Journal of Neuroscience*, 42(27). [doi:10.1523/JNEUROSCI.2203-21.2022](https://doi.org/10.1523/JNEUROSCI.2203-21.2022)

Humans and macaques steer with a joystick through a virtual ground plane to a remembered "firefly"
location, using only the optic flow their own movement generates. To test whether they integrate that
flow over time (rather than replaying a learned motor plan), the paper uses three manipulations, and
this repository holds the analysis for each:

| Module | Manipulation | What the code does |
|---|---|---|
| `perturbation/` | Unpredictable optic-flow perturbations that push the subject off course | The two core methods: `gen_sim_ptb` simulates the "ignored the perturbation" counterfactual by adding the perturbation velocity profile to matched unperturbed trajectories and integrating to a stopping point (the control against which real endpoints are compared by ROC/AUC); `get_kernel` estimates the velocity response kernel to the perturbation by regressing speed traces on an amplitude-scaled boxcar basis |
| `gain/` | Joystick gain (1, 1.5, 2) changes the consequence of every action | Pools sessions by gain (`con_gain`); travel-time CDFs per gain (`time_cdf_fun`, `time_cdf`); log-log regression of travel time on distance and speed (`LOGregress_fun`) |
| `density/` | Optic-flow density changes how much evidence the flow carries | Pools sessions by ground-plane density and splits endpoints per density level (`con_dens`) |
| `utils/` | | Euclidean endpoint error, cell concatenation |

## How it fits with the lab pipeline

Raw sessions are imported and pre-processed by the Angelaki lab's
[`firefly-monkey`](https://github.com/kaushik-l/firefly-monkey) pipeline (Kaushik Lakshminarasimhan),
which produces one `behaviours` object per session (trials with `events`, `prs`, `logical`,
`continuous`, and `stats` with final positions and trial-type indices). Everything here consumes those
objects; nothing from that pipeline is duplicated in this repository. `ComputeROCFirefly` (ROC
analysis) also comes from `firefly-monkey` and must be on the path.

## Running

MATLAB R2020a or later, Statistics and Machine Learning Toolbox (`regress`), Image Processing Toolbox
(`medfilt3`), and `firefly-monkey` on the path. Point the scripts at your data:

```matlab
setenv('FIREFLY_DATA_ROOT', '/path/to/sessions')   % folders m44/, m51/, m53/ with session .mat files
addpath(genpath(pwd))
time_cdf                                            % gain: travel-time CDFs per gain, all animals
m = con_dens(fullfile(getenv('FIREFLY_DATA_ROOT'),'m51'));   % density: pooled endpoints per density
[kernel, ts] = get_kernel(Yresp, Amplitudes, Yphant);        % perturbation: response kernel
```

Behavioural data are not distributed here; see the data statement in the paper.

## Licence

MIT (see `LICENSE`).
