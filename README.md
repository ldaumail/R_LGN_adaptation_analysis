# R_LGN_adaptation_analysis
### Statistical analyses of the rapid visual adaptation manuscript were performed on R
## This repository hosts the following analyses:
## -Adaptation analysis:
## *LMMs fit of the peaks and troughs (individual_channels_analysis_origpeaks_KRpvalcorr.R)
             -data location:- "C:/Users/daumail/OneDrive -        Vanderbilt/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/*.mat"
## *LMMs fit on peaks in both mono and binocular conditions:(lmer_peaks_binocular_adaptation.R)
             -data location:- "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/all_units/all_orig_bs_zscore_trials_05022021_mono_bino.mat"
## *Some t-tests(dependent_t_test_normality_test.R)
             -data location:- "/Users/loicdaumail/Documents/Research_MacBook/Maier_lab/adaptation_LGN_project/R/all_raw_mean_data_peaks.mat"
## *Bootstrap analysis (bootstrap_anal.R)
             -data location:- "/Users/loicdaumail/Documents/Research_MacBook/Maier_lab/adaptation_LGN_project/R/all_raw_mean_data_peaks.mat"
## *Effect size (cohens_d_plots.R)
             -data location:-"/Users/loicdaumail/Documents/Research_MacBook/Maier_lab/adaptation_LGN_project/R/all_raw_mean_data_peaks.mat"
## *plots(violin_plots_edited.R)
             -data location:-"/Users/loicdaumail/Documents/Research_MacBook/Maier_lab/adaptation_LGN_project/R/all_norm_mean_data_peaks.mat"
## -Binocular interaction analysis:
## *LMs fit of linear trends/anovas and effect sizes =proportion of variance explained regarding each type of trend (contrast_linearTrend_and_monovsBino_interaction.R)
              -data location:- "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data/single_units/binocular_adaptation/all_units/all_orig_bs_zscore_trials.mat"
## *Other attempts to fit better models to the data to test for interaction (all_Trends_and_monovsBino_interaction.R,contrast_Pk1vsPk4_and_Mono_vs_Bino_interaction.R)
