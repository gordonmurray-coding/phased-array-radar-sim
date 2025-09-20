function run_preset(presetName)
%RUN_PRESET Convenience wrapper to run the radar sim with a named preset.
%   run_preset('high_res_range')

if nargin < 1, presetName = 'default'; end
Radar_theoretical_max_accuracy_formulas_waveforms('preset', presetName);
end
