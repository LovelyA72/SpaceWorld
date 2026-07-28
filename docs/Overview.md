# SpaceWorld4U Docs

> Hey there! This is **Kashouryo** (何小綾), Lovely's lab assistant, and you can call me Shouryo, or just **KSR**. I’m certainly not a singer myself (But I am working on it... hopefully!), but I’ll be guiding you through this document. It's sooooo nice to meet you! :ksr_suki:

SpaceWorld4U(SW4U or SpaceWorld) is a OpenUTAU/UTAU resampler.It takes a short WAV or AIFF recording, changes its pitch and length, and writes a mono WAV file.

SW4U's signal processing is based on the WORLD vocoder. This repository packages the vocoder and the resampler in one C++ executable. No need to mess with external libraries. (You're welcome)

> KSR here! If the input doesn't have a `.frq` file yet, SW4U generates one automatically before running the rest of the analysis.

## How it works

Most SW4U behavior is in `main.cpp`. A normal run follows these steps:

1. Look beside the input audio file for cached WORLD analysis data.
2. If the cache is missing, invalid, or the `R` flag is present, read the input audio and its UTAU `.frq` frequency file. If neither supported `.frq` name exists, generate one first.
3. Analyze the voice:
   - DIO estimates the fundamental frequency (F0).
   - data from the UTAU `.frq` file repairs large DIO errors.
   - StoneMask refines the F0 values.
   - CheapTrick calculates the spectral envelope, which represents vocal tone.
   - D4C calculates aperiodicity, which represents noise and breathiness.
4. Map input frames to output frames. This is where offset, fixed region, cutoff, velocity, requested length, note, modulation, pitch bend, and the looping option are applied.
5. Optionally change the spectral envelope with the `g` flag, then synthesize a new waveform with WORLD.
6. Apply the `B` and `O` post-processing flags, normalize/scale the waveform, and write the output WAV file.

The analysis frame period is fixed in `main.cpp` at 256 samples at 44.1 kHz, which is about 5.805 milliseconds. Several calculations and cache formats assume this value.

> KSR here! That 5.805 ms frame period is baked into timing and cache behavior. If you change it, please treat the cache format and every frame-to-time calculation as part of the same change!

## Repository map

### Root files

| Path | Purpose |
| --- | --- |
| `README.md` | Very short project introduction and CI badge. |
| `SpaceWorld.sln` | Visual Studio solution. |
| `SpaceWorld/SpaceWorld.vcxproj` | Source list and Debug/Release, x86/x64 build settings. |
| `.github/workflows/msbuild.yml` | Builds the Release configuration on Windows after each push and uploads the x64 executable. It does not run tests. |
| `LICENSE` | GPL-3.0 license text. |

### Application and audio code

| Path | Purpose |
| --- | --- |
| `SpaceWorld/src/main.cpp` | CLI parsing, cache handling, analysis orchestration, time/pitch transformation, effects, synthesis, and WAV writing. Start here when tracing program behavior. |
| `SpaceWorld/src/audio_io.cpp/.h` | Active WAV/AIFF reader. Stereo input is averaged to mono. |
| `SpaceWorld/src/wavread.cpp/.h` | Older WAV reader. It is compiled but is not used by the current main path. |
| `SpaceWorld/src/getWorldValues77.cpp/.h` | Frame, sample, time, frequency, and FFT-size helper functions. |
| `SpaceWorld/src/matlabfunctions.cpp` | Math and interpolation helpers originally modeled on MATLAB functions. |
| `SpaceWorld/src/matlabmyfunctions.cpp/.h` | More interpolation, windowing, sorting, and numeric helpers. |
| `SpaceWorld/src/world.h` | Older shared WORLD declarations and project-wide constants such as F0 limits. |

### WORLD vocoder code

| Path | Purpose |
| --- | --- |
| `SpaceWorld/src/dio.cpp` | DIO F0 estimation. |
| `SpaceWorld/src/stonemask.cpp` | StoneMask F0 refinement. |
| `SpaceWorld/src/cheaptrick.cpp` | CheapTrick spectral-envelope analysis. |
| `SpaceWorld/src/d4c.cpp` | D4C aperiodicity analysis. |
| `SpaceWorld/src/synthesis.cpp` | Current WORLD waveform synthesis. |
| `SpaceWorld/src/common.cpp` | Shared spectrum, smoothing, window, and FFT setup code. |
| `SpaceWorld/src/fft.cpp` | Bundled Ooura FFT implementation and WORLD-compatible wrappers. |
| `SpaceWorld/src/world/*.h` | Public declarations and constants used by the newer WORLD source files. |
| `SpaceWorld/src/star.cpp` | Older STAR analysis code. The current runtime uses its FFT-size helper but uses CheapTrick for analysis. |
| `SpaceWorld/src/platinum.cpp` | Older PLATINUM residual analysis code. It is compiled but is not called by the current main path. |

`dio(old).cpp` and `synthesis (old).cpp` are kept as reference snapshots. They are deliberately not listed in the Visual Studio project and are not compiled.

## Building

Requirements:

- Windows
- Visual Studio 2022 with the **Desktop development with C++** workload, or equivalent MSBuild v143 tools
- A Windows 10 SDK

Open `SpaceWorld.sln` in Visual Studio and build the wanted configuration, or run this from a Visual Studio Developer PowerShell:

```powershell
msbuild /m /p:Configuration=Release /p:Platform=x64 SpaceWorld.sln
```

The x64 post-build step renames the executable to:

```text
x64/Release/SpaceWorld_win64.exe
```

Debug, Release, x86, and x64 configurations exist. The source has direct Windows dependencies, including `windows.h`, `timeGetTime`, and `winmm.lib`, so it is not currently a portable Linux/macOS project.

There are no NuGet packages or other external dependencies in the repository. The GitHub workflow still runs `nuget restore`, but there is no package manifest to restore.

## Command-line interface

The program uses the positional interface expected by a UTAU resampler:

```text
SpaceWorld_win64.exe input output note velocity flags offset length fixed cutoff volume modulation tempo pitchbend
```

Only the first four values after the executable are required by the argument count check. Later values are optional in C++ but are positional, so callers must provide placeholders for earlier values when they want to set a later one.

> KSR here! These arguments are positional, so don’t skip an earlier one just because you only want a later option! Give it a placeholder, or every argument after it will slide into the wrong slot. Total chaos~

| Position | Name | Meaning |
| ---: | --- | --- |
| 1 | `input` | Input WAV or AIFF path. |
| 2 | `output` | Output WAV path. |
| 3 | `note` | Target note such as `C4`, `F#4`, or `Bb3`. |
| 4 | `velocity` | UTAU velocity value. `100` is the neutral value. |
| 5 | `flags` | One combined string containing zero or more flags described below. |
| 6 | `offset` | Milliseconds skipped at the start of the input. |
| 7 | `length` | Requested output length in milliseconds. |
| 8 | `fixed` | Leading region, in milliseconds, handled separately from the stretchable region. |
| 9 | `cutoff` | Positive: milliseconds removed from the end. Negative: usable length measured forward from the offset. |
| 10 | `volume` | Output volume percent; defaults to 100. |
| 11 | `modulation` | Percent of the source F0 movement kept in the result; defaults to 100. |
| 12 | `tempo` | UTAU tempo token, normally like `!120`. The code skips its first character. |
| 13 | `pitchbend` | UTAU's compact, base-64-like pitch-bend string. Values are cents and are interpolated over time. |

Supported flags in position 5:

| Flag | Range/default | Effect |
| --- | --- | --- |
| `R` | off | Ignore cached analysis and analyze the input again. |
| `tN` | `0` | Shift the target note by `N` cents. |
| `e` | off | Loop the stretchable region back and forth instead of only stretching it. |
| `gN` | -100 to 100 | Move the spectral envelope to change the apparent vocal size/timbre. |
| `BN` | 0 to 100; default 50 | Change breath/noise character. Above 50 adds noise; below 50 applies a comb filter. |
| `ON` | -100 to 100; default 0 | Change vocal opening/strength with spectral post-processing. |
| `PN` | 0 to 100; default 86 | Control peak compensation during output normalization. |
| `G` | off | Parsed by `main.cpp`, but currently has no effect. |

Example:

```powershell
.\x64\Release\SpaceWorld_win64.exe input.wav output.wav C4 100 "g-10B60" 0 500 100 0 100 100
```

If compatible analysis caches and a UTAU frequency file are both missing,
SpaceWorld generates the frequency file before rendering.

> KSR here, again! Capital G is recognized but doesn’t actually do anything yet. (You don't want me, such a cute lady, to growl! Don't ya?) Lowercase gN is the working timbre flag, so mind the capitalization!

## Input, output, and sidecar files

### Audio

The active reader accepts:

- WAV: mono or stereo; 8-, 16-, or 24-bit PCM, or 32-bit IEEE float
- AIFF: mono or stereo; signed 8-, 16-, or 24-bit PCM

Stereo is mixed down by averaging the left and right channels.

The writer always stores mono 16-bit sample values. However, it currently puts the input bit depth into parts of the WAV header. In practice, use 16-bit input when testing normal resampler behavior until this header bug is fixed.

> KSR here! Please use 16-bit input for now, okay? SW4U always writes 16-bit samples, but the header still copies the input bit depth. Using 24-bit or float input may create a weird WAV file! :ksr_nooo:

### UTAU frequency file

The reader accepts version `FREQ0003` and looks for these names, in order:

1. the full input name plus `.frq`, for example `voice.wav.frq`
2. the input extension changed to an underscore plus `.frq`, for example `voice_wav.frq`

If neither file is found, SpaceWorld analyzes the audio and writes the second
form, such as `voice_wav.frq`. An existing but invalid `.frq` is not
overwritten. If generation or subsequent loading fails, the program exits
with code `-60`.

### WORLD analysis cache

After analysis, the program writes three binary files beside the input. The input extension is replaced:

| Extension | Contents |
| --- | --- |
| `.dio` | Time positions and refined F0 values. |
| `.ctspec` | Compressed CheapTrick spectral envelopes. |
| `.d4c` | Compressed D4C aperiodicity values. |

For `voice.wav`, these are `voice.dio`, `voice.ctspec`, and `voice.d4c`.

All three must load successfully for the fast cached path. The cache checks stored dimensions, sample rate, and input sample count, but it does not compare file timestamps or a content hash. Use the `R` flag after changing audio while keeping the same length and sample rate.

> KSR here! Changed the source audio but kept the same sample rate and length? Use the R flag! The cache can’t tell that the actual audio changed, so it may happily reuse stale analysis. Sneaky, right?

These generated files are not covered by a special rule in `.gitignore`. Avoid committing voice data and analysis caches unless they are intentional test fixtures.

## Testing and validation

> KSR here! We don’t have automated tests yet, so a successful build isn’t enough. Please listen to the result and check its duration, pitch, silence, clipping, and stretch boundaries before declaring victory!

There is currently no automated test suite. CI only confirms that the Release solution builds on `windows-latest`.

For a signal-processing change, a useful manual check is:

1. Build Release x64.
2. Run once with a known 16-bit voice sample and no `.frq`.
3. Confirm that `.frq`, `.dio`, `.ctspec`, and `.d4c` are created.
4. Run the same command again and confirm that the cache is read.
5. Run with `R` and confirm that analysis is rebuilt.
6. Inspect the WAV duration, header, clipping, silence, pitch, and audible transitions around the fixed/stretch boundary.
7. Repeat with the specific flags changed by the patch.

When adding automated tests, small tests for note parsing, pitch-bend decoding, time mapping, cache headers, and WAV headers would give the highest early value. A short licensed audio fixture could then support one end-to-end test.

## Important maintenance notes

- `main.cpp` is almost 2,000 lines and owns most of the application logic. Small, pure helper functions are the safest first extraction targets.
- Memory management mixes `malloc/free` and `new[]`; some paths do not use a matching deallocator. Treat memory-lifetime changes carefully and use runtime checking when possible.
- Many errors return process code `0`, while a failed `.frq` generation or load calls `exit(-60)`. Do not assume that `0` always means success until error handling is cleaned up.
- File names are copied into fixed-size C buffers. Long paths and unchecked writes are a known risk.
- Some comments are Japanese, and some older comments have damaged character encoding. Check the code itself before relying on an unclear comment.
- New `.cpp` or `.h` files must be added to `SpaceWorld.vcxproj`; simply placing them under `src` does not make Visual Studio compile them.
- The current and old implementations have similar symbols. Do not add `dio(old).cpp` or `synthesis (old).cpp` to the build unless the duplicate definitions are intentionally resolved.
- Cache formats are raw native C/C++ values with no general versioning layer. Changing types, dimensions, frame timing, or compression requires a cache compatibility plan.
- Signal-processing arrays use `0.0` F0 for unvoiced frames. Preserve that convention when changing interpolation or pitch logic.

For most behavior changes, begin in `main.cpp`, then follow the called function into the matching WORLD module. For algorithm changes, keep CLI, cache, and file-format work separate from the DSP change so that regressions are easier to locate.
