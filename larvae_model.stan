// Multivariate mixed model for larval traits
// Traits: Trunk length, Tail length (continuous MVN),
//         Hatching probability, Settling|Hatching probability (Bernoulli/logit)
//
// Random effects: sire, dam, sire-x-dam interaction
//   each drawn from a 4-variate normal with estimated covariance matrix
// Fixed effects: block means (one per block per trait)
//
// Changes from previous version:
//   1. SD priors: half-normal(0, sqrt(1e3)) replacing half-uniform(0, sqrt(1e3)).
//      Equivalent upper mass but curvature near zero prevents funnel geometry.
//   2. LKJ(2) replacing LKJ(1) on all correlation matrices.
//      Regularises away from near-singular correlations.
//   3. Trunk/tail likelihood uses multi_normal_cholesky() on a pre-built
//      mean matrix, avoiding per-observation Cholesky decompositions.

data {
  // Dimensions
  int<lower=1> n_blocks;
  int<lower=1> n_sires;
  int<lower=1> n_dams;
  int<lower=1> n_int;

  // Trunk/tail data
  int<lower=1> n_larvae;
  array[n_larvae] int<lower=1, upper=n_blocks> tt_block;
  array[n_larvae] int<lower=1, upper=n_sires>  tt_sire;
  array[n_larvae] int<lower=1, upper=n_dams>   tt_dam;
  array[n_larvae] int<lower=1, upper=n_int>    tt_int;
  matrix[n_larvae, 2] length_obs;   // columns: Trunk, Tail

  // Hatching data
  int<lower=1> n_hatch;
  array[n_hatch] int<lower=1, upper=n_blocks> h_block;
  array[n_hatch] int<lower=1, upper=n_sires>  h_sire;
  array[n_hatch] int<lower=1, upper=n_dams>   h_dam;
  array[n_hatch] int<lower=1, upper=n_int>    h_int;
  array[n_hatch] int<lower=0, upper=1> hatch;

  // Settling|Hatching data
  int<lower=1> n_settle_hatch;
  array[n_settle_hatch] int<lower=1, upper=n_blocks> sh_block;
  array[n_settle_hatch] int<lower=1, upper=n_sires>  sh_sire;
  array[n_settle_hatch] int<lower=1, upper=n_dams>   sh_dam;
  array[n_settle_hatch] int<lower=1, upper=n_int>    sh_int;
  array[n_settle_hatch] int<lower=0, upper=1> settle_hatch;

  // Per-trait prior bounds for block means [2 x 4]: row 1 = lower, row 2 = upper
  // Columns: Trunk, Tail, Hatch (logit scale), Settle|Hatch (logit scale)
  matrix[2, 4] block_mean_range;
}

parameters {
  // ---- block means: unit-interval raw, rescaled in transformed parameters ----
  matrix<lower=0, upper=1>[n_blocks, 4] block_mean_raw;

  // ---- random effect means ----
  vector[4] sire_mean;
  vector[4] dam_mean;
  vector[4] int_mean;

  // ---- standard deviations ----
  // half-normal(<lower=0> declaration + normal(0,s) statement = half-normal prior)
  // Scale sqrt(1e3) ≈ 31.6 matches the original var ~ Uniform(0, 1e3) upper bound
  // while providing curvature near zero that prevents funnel geometry.
  vector<lower=0>[4] sire_sd;
  vector<lower=0>[4] dam_sd;
  vector<lower=0>[4] int_sd;
  vector<lower=0>[2] resid_sd;

  // ---- Cholesky factors of correlation matrices ----
  // LKJ(2): regularises away from near-singular correlations.
  cholesky_factor_corr[4] L_sire;
  cholesky_factor_corr[4] L_dam;
  cholesky_factor_corr[4] L_int;
  cholesky_factor_corr[2] L_resid;

  // ---- non-centred random effect z-scores ----
  matrix[n_sires, 4] z_sire;
  matrix[n_dams,  4] z_dam;
  matrix[n_int,   4] z_int;
}

transformed parameters {
  // Rescale block means from [0,1] into their prior support per trait
  matrix[n_blocks, 4] block_mean;
  for (t in 1:4)
    block_mean[, t] = block_mean_range[1, t] +
      block_mean_raw[, t] * (block_mean_range[2, t] - block_mean_range[1, t]);

  // Cholesky factors of covariance matrices (used directly in likelihood)
  matrix[4, 4] L_sire_vcov  = diag_pre_multiply(sire_sd,  L_sire);
  matrix[4, 4] L_dam_vcov   = diag_pre_multiply(dam_sd,   L_dam);
  matrix[4, 4] L_int_vcov   = diag_pre_multiply(int_sd,   L_int);
  matrix[2, 2] L_resid_vcov = diag_pre_multiply(resid_sd, L_resid);

  // Full covariance matrices (for posterior output / QG calculations)
  matrix[4, 4] sire_vcov  = L_sire_vcov  * L_sire_vcov';
  matrix[4, 4] dam_vcov   = L_dam_vcov   * L_dam_vcov';
  matrix[4, 4] int_vcov   = L_int_vcov   * L_int_vcov';
  matrix[2, 2] resid_vcov = L_resid_vcov * L_resid_vcov';

  // Realised random effects (non-centred parameterisation)
  matrix[n_sires, 4] sire_eff;
  matrix[n_dams,  4] dam_eff;
  matrix[n_int,   4] int_eff;

  for (s in 1:n_sires)
    sire_eff[s] = sire_mean' + (L_sire_vcov * z_sire[s]')';
  for (d in 1:n_dams)
    dam_eff[d]  = dam_mean'  + (L_dam_vcov  * z_dam[d]')';
  for (i in 1:n_int)
    int_eff[i]  = int_mean'  + (L_int_vcov  * z_int[i]')';
}

model {
  // ---- priors ----

  // block_mean_raw ~ Uniform(0,1) implicit from <lower=0,upper=1> declaration.

  // Effect means: precision 1e-3 matches JAGS dnorm(0, 1e-3)
  sire_mean ~ normal(0, sqrt(1e3));
  dam_mean  ~ normal(0, sqrt(1e3));
  int_mean  ~ normal(0, sqrt(1e3));

  // SDs: half-normal, scale = sqrt(1e3) ≈ 31.6
  // Equivalent upper support to original Uniform(0, sqrt(1e3)) but with
  // curvature near zero; avoids funnel geometry in the hierarchical posterior.
  sire_sd  ~ normal(0, sqrt(1e3));
  dam_sd   ~ normal(0, sqrt(1e3));
  int_sd   ~ normal(0, sqrt(1e3));
  resid_sd ~ normal(0, sqrt(1e3));

  // Correlation matrices: LKJ(2) concentrates away from near-singular matrices.
  L_sire  ~ lkj_corr_cholesky(2);
  L_dam   ~ lkj_corr_cholesky(2);
  L_int   ~ lkj_corr_cholesky(2);
  L_resid ~ lkj_corr_cholesky(2);

  // Non-centred z-scores
  to_vector(z_sire) ~ std_normal();
  to_vector(z_dam)  ~ std_normal();
  to_vector(z_int)  ~ std_normal();

  // ---- likelihood: trunk/tail (bivariate normal) ----
  // Vectorised over larvae: build full mean matrix first, then call
  // multi_normal_cholesky() once per observation using the already-computed
  // L_resid_vcov. Avoids n_larvae redundant Cholesky decompositions.
  {
    matrix[n_larvae, 2] mu;
    for (l in 1:n_larvae)
      mu[l] = block_mean[tt_block[l], 1:2] +
               sire_eff[tt_sire[l], 1:2]   +
               dam_eff[tt_dam[l],   1:2]   +
               int_eff[tt_int[l],   1:2];
    for (l in 1:n_larvae)
      length_obs[l] ~ multi_normal_cholesky(mu[l]', L_resid_vcov);
  }

  // ---- likelihood: hatching (Bernoulli logistic) ----
  {
    vector[n_hatch] eta_h;
    for (h in 1:n_hatch)
      eta_h[h] = block_mean[h_block[h], 3] +
                 sire_eff[h_sire[h], 3]    +
                 dam_eff[h_dam[h],   3]    +
                 int_eff[h_int[h],   3];
    hatch ~ bernoulli_logit(eta_h);
  }

  // ---- likelihood: settling|hatching (Bernoulli logistic) ----
  {
    vector[n_settle_hatch] eta_sh;
    for (sh in 1:n_settle_hatch)
      eta_sh[sh] = block_mean[sh_block[sh], 4] +
                   sire_eff[sh_sire[sh], 4]    +
                   dam_eff[sh_dam[sh],   4]    +
                   int_eff[sh_int[sh],   4];
    settle_hatch ~ bernoulli_logit(eta_sh);
  }
}

generated quantities {
  // ---- overall mean and variance (Trunk/Tail) for evolvability calculations ----
  vector[2] overall_mean;
  vector[2] overall_var;
  {
    vector[2] sum_y  = rep_vector(0.0, 2);
    vector[2] sum_y2 = rep_vector(0.0, 2);
    for (l in 1:n_larvae) {
      sum_y  += length_obs[l]';
      sum_y2 += square(length_obs[l]');
    }
    overall_mean = sum_y  / n_larvae;
    overall_var  = sum_y2 / n_larvae - square(overall_mean);
  }

  // ---- posterior predictive draws ----
  matrix[n_larvae, 2] length_ppd;
  array[n_hatch] int hatch_ppd;
  array[n_settle_hatch] int settle_hatch_ppd;

  {
    for (l in 1:n_larvae) {
      vector[2] mu = block_mean[tt_block[l], 1:2]' +
                     sire_eff[tt_sire[l], 1:2]'   +
                     dam_eff[tt_dam[l],   1:2]'   +
                     int_eff[tt_int[l],   1:2]';
      length_ppd[l] = to_row_vector(
        multi_normal_cholesky_rng(mu, L_resid_vcov)
      );
    }
  }

  for (h in 1:n_hatch) {
    real eta = block_mean[h_block[h], 3] +
               sire_eff[h_sire[h], 3]   +
               dam_eff[h_dam[h],   3]   +
               int_eff[h_int[h],   3];
    hatch_ppd[h] = bernoulli_logit_rng(eta);
  }

  for (sh in 1:n_settle_hatch) {
    real eta = block_mean[sh_block[sh], 4] +
               sire_eff[sh_sire[sh], 4]   +
               dam_eff[sh_dam[sh],   4]   +
               int_eff[sh_int[sh],   4];
    settle_hatch_ppd[sh] = bernoulli_logit_rng(eta);
  }
}
