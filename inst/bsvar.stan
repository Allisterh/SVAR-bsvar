functions {
  real sged_lpdf(vector x, int N, real lambda, real log_p, int fix_moments) {
    real p = exp(log_p);
    real mu = 0;
    real sigma = 1;
    real m = 0;
    real v = 1;
    real f1;
    real f2;
    vector[N] xm;
    vector[N] abs_xm;
    vector[N] ret;
    if (fix_moments > 1) {
      v = pow( (pi() * (1 + 3 * lambda^2) * exp(lgamma(3/p)) - 16^(1/p) * lambda^2 * (exp(lgamma(0.5 + (1/p))))^2 * exp(lgamma(1/p))) / (pi() * exp(lgamma(1/p))), (-0.5));
    }
    if (fix_moments > 0) {
      m = (2^(2/p) * v * sigma * lambda * exp(lgamma(0.5 + (1/p)))) / sqrt(pi());
    }
    f1 = log(p) - log(2) - log(v) - log(sigma) - lgamma(1/p);
    xm = x + m - mu;
    abs_xm = fabs(xm);
    for(i in 1:N) {
      if ((xm[i]) < 0) {
        f2 = (v * sigma * (1 - lambda));
      } else {
        f2 = (v * sigma * (1 + lambda));
      }
      ret[i] = f1 - (abs_xm[i] / f2)^p;
    }
    return sum(ret);
  }
  real sgt_lpdf(vector x, int N, real lambda, real log_p, real log_q, int fix_moments) {
    real p = exp(log_p);
    real q = exp(log_q);
    real mu = 0;
    real sigma = 1;
    real m = 0;
    real v = 1;
    real ipq = (1/p + q);
    real f1;
    real f2_base;
    real f2_neg;
    real f2_pos;
    vector[N] xm;
    vector[N] abs_xm;
    vector[N] ret;
    if(log_q == positive_infinity()) {
      return(sged_lpdf(x | N, lambda, log_p, fix_moments));
    }
    if (fix_moments > 1) {
      v = pow(q, - (1/p)) * pow((3 * pow(lambda, 2) + 1) * (exp(lbeta((3/p), (q - (2/p)))) / exp(lbeta((1/p), q))) - 4 * pow(lambda, 2) * pow((exp(lbeta((2/p), (q - (1/p)))) / exp(lbeta((1/p), q))), 2), -(0.5));
    }
    if (fix_moments > 0) {
      m = (2 * v * sigma * lambda * pow(q, (1/p)) * exp(lbeta((2/p), (q-(1/p))))) / exp(lbeta((1/p), q));
    }
    f1 = log_p - log(2) - log(v) - log(sigma) - (log_q/p) - lbeta((1/p), q);
    f2_base = q * pow((v * sigma), p);
    f2_neg = (f2_base * pow((lambda*(-1) + 1), p));
    f2_pos = (f2_base * pow((lambda + 1), p));
    xm = x + m - mu;
    abs_xm = fabs(xm);
    for(i in 1:N) {
      real log1p_arg;
      if ((xm[i]) < 0) {
        log1p_arg = pow(abs_xm[i], p) / f2_neg;
      } else {
        log1p_arg = pow(abs_xm[i], p) / f2_pos;
      }
      ret[i] = f1 - ipq * log1p(log1p_arg);
    }
    return sum(ret);
  }
  real likelihood_partial_sum(array[] int col_indices, int start, int end, int N, matrix shocks, matrix gamma, int fix_moments) {
    real out = 0;
    for(i in start:end) {
      out += sgt_lpdf(col(shocks, i) | N, gamma[i, 1], gamma[i, 2], gamma[i, 3], fix_moments);
    }
    return out;
  }
  real likelihood_vol_partial_sum(array[] int col_indices, int start, int end, int N, matrix shocks_init, matrix gamma, int fix_moments,
                                  array[] vector garch_param, vector garch_v_init, int garch_eta_form) {
    matrix[rows(shocks_init), cols(shocks_init)] shocks = shocks_init;
    vector[rows(garch_v_init)] garch_v = garch_v_init;
    real out = 0;
    for(j in start:end) {
      int vol_break_count = 1;
      for(i in 2:N) {
        if(garch_eta_form == 0) garch_v[j] = sqrt(garch_param[j][1] + garch_param[j][2] * garch_v[j]^2 + garch_param[j][3] * (shocks[i-1,j] * garch_v[j])^2);
        if(garch_eta_form == 1) garch_v[j] = sqrt(garch_param[j][1] + garch_param[j][2] * garch_v[j]^2 + garch_param[j][3] * shocks[i-1,j]^2);
        shocks[i,j] = shocks[i,j] / garch_v[j];
        out += -log(fabs(garch_v[j]));
      }
      out += sgt_lpdf(sub_col(shocks, 2, j, N - 1) | N - 1, gamma[j, 1], gamma[j, 2], gamma[j, 3], fix_moments);
    }
    return out;
  }
  real minnesota_partial_sum(array[] int col_indices, int start, int end, int lags, int M, matrix A, array[] real minnesota_means, array[] real hypers, vector reduced_scales, array[] int no_cross_shrinkage_vec, array[,] int A_zero_res) {
    real out = 0;
    for(i in start:end) {
        for(h in 1:lags) {
          for(j in 1:M) {
            if(A_zero_res[(j-1)*M*lags + (h-1)*M + i, 1] != 1) {
              if (i == j) {
                if (h == 1) {
                  out += normal_lpdf(A[(h-1)*M + i, j] | minnesota_means[i], hypers[1]);
                } else {
                  out += normal_lpdf(A[(h-1)*M + i, j] | 0, hypers[1] / (h^hypers[2]));
                }
              } else {
                if(no_cross_shrinkage_vec[i] + no_cross_shrinkage_vec[j] == 2) {
                  out += normal_lpdf(A[(h-1)*M + i, j] | 0, (hypers[1] * reduced_scales[j]) / ((h^hypers[2]) * reduced_scales[i]) );
                } else {
                  out += normal_lpdf(A[(h-1)*M + i, j] | 0, (hypers[1] * hypers[3] * reduced_scales[j]) / ((h^hypers[2]) * reduced_scales[i]) );
                }
              }
            }
          }
        }
      }
    return out;
  }
  matrix build_yx(matrix y_raw, int lags, vector missing_data, array[,] int missing_data_location, int missing_data_count, int N, int M) {
    matrix[N + lags, M] y_new = y_raw;
    matrix[N + lags, M*lags] xx_raw;
    matrix[N + lags, M + M*lags] yx_raw;
    if(missing_data_count > 0) {
      for(i in 1:missing_data_count) {
        y_new[missing_data_location[i, 1], missing_data_location[i, 2]] = y_new[(missing_data_location[i, 1] - 1), missing_data_location[i, 2]] + missing_data[i];
      }
    }
    for(i in 1:lags) {
      for(j in 1:M) {
        int xx_col = (i - 1)*M + j;
        for(xx_row in (i + 1):(N + lags)) {
          xx_raw[xx_row, xx_col] = y_new[(xx_row - i), j];
        }
      }
    }
    yx_raw[,1:M] = y_new;
    yx_raw[,(M + 1):(M + M*lags)] = xx_raw;
    return yx_raw;
  }
  matrix rebuild_y_raw(matrix y_raw, array[] int quarterly_binary, vector quarterly_sums, array[] vector monthly_raw, vector data_transformation_mean, vector data_transformation_scale, int M, int NQ) {
    matrix[rows(y_raw), M] y_new = y_raw;
    int count = 1;
    for(j in 1:M) {
      if(quarterly_binary[j] == 1) {
        for(i in 1:NQ) {
          vector[3] sub = log(quarterly_sums[count] * monthly_raw[count]); // Log-levels
          sub = (sub - data_transformation_mean[j]) / data_transformation_scale[j]; // Data scaling
          y_new[(i * 3 - 2):(i * 3), j] = sub;
          count += 1;
        }
      }
    }
    return y_new;
  }
  matrix subtract_measurement_errors(matrix y_raw, array[] int measurement_errors_binary, vector measurement_errors, int M, int N_plus_lags) {
    matrix[rows(y_raw), M] y_new = y_raw;
    int me_count = 0;
    for(j in 1:M) {
      if(measurement_errors_binary[j] == 1) {
        for(i in 1:N_plus_lags) {
          if(y_raw[i,j] != positive_infinity()) {
            me_count += 1;
            y_new[i,j] = y_new[i,j] - measurement_errors[me_count];
          }
        }
      }
    }
    return y_new;
  }
}
data {
  int<lower=0> N;
  int<lower=1> M;
  int<lower=0> lags;
  matrix[N, M] y;
  matrix[N, M*lags] x;
  matrix[(N == 0 ? N : N + lags), M] y_raw;
  int<lower=0, upper=1> include_constant;
  array[2] real constant_prior;
  int<lower=0, upper=1> B_inverse;
  int<lower=0> pq_len;
  array[2] real<lower=0> p_q_mins;
  array[2] real lambda_prior;
  array[2] real p_prior;
  array[2] real q_prior;
  array[2] real p_q_prior_shift;
  array[3] int sgt_len_vec;
  array[3] real sgt_fixed;
  int<lower=0, upper=M*M*lags> A_zero_res_num;
  array[M*M*lags, 2] int A_zero_res;
  int<lower=0, upper=1> include_minnesota;
  array[3] int hyper_len_vec;
  array[3] real hyper_fixed;
  array[2] real hyper_shrinkage_prior;
  array[2] real hyper_lags_prior;
  array[2] real hyper_ownlags_prior;
  array[2] real<lower=0> hyper_shrinkage_lim;
  int<lower=0, upper=1> include_soc;
  int<lower=0, upper=1> include_dio;
  int<lower=0, upper=1> soc_free;
  int<lower=0, upper=1> dio_free;
  real soc;
  real dio;
  array[2] real soc_prior;
  array[2] real dio_prior;
  array[2] int soc_dio_dep;
  array[M] real minnesota_means;
  int<lower=0, upper=2> fix_moments;
  int<lower=0, upper=M*M> B_zero_res_num;
  int<lower=0, upper=M*M> B_pos_res_num;
  int<lower=0, upper=M*M> B_neg_res_num;
  array[M*M, 2] int B_zero_res;
  array[M*M, 2] int B_pos_res;
  array[M*M, 2] int B_neg_res;
  int<lower=0, upper=1> include_B_prior;
  vector[M*M] B_prior_mean;
  matrix[M*M, M*M] B_prior_cov;
  int<lower=0, upper=1> B_prior_cov_diagonal;
  int<lower=0, upper=1> include_garch;
  int<lower=0, upper=1> garch_dependence;
  vector<lower=0>[include_garch * (garch_dependence == 1 ? M + 2 : 3)] garch_prior;
  int<lower=0, upper=1> include_garch_groups;
  int<lower=0> garch_group_num;
  matrix[include_garch_groups * M, garch_group_num] garch_group_mat;
  int<lower=0, upper=1> garch_eta_form;
  int<lower=0, upper=1> include_vol_breaks;
  int<lower=0> vol_breaks_num;
  array[N] int<lower=0, upper=1> vol_breaks;
  vector<lower=0>[vol_breaks_num + 1] vol_breaks_prior;
  int<lower=0, upper=1> include_hyper_B;
  vector<lower=0>[2] hyper_B_prior;
  array[M * M] int<lower=0, upper=1> prior_elast_binary;
  int<lower=0> prior_elast_num;
  array[prior_elast_num, 2] int prior_elast_ind;
  matrix[prior_elast_num, 3] prior_elast_par;
  int<lower=0, upper=1> include_missing_data;
  int<lower=0> missing_data_count;
  array[missing_data_count * include_missing_data, 2] int missing_data_location;
  int<lower=0> number_of_quarterly;
  array[M] int<lower=0, upper=1> quarterly_binary;
  vector[number_of_quarterly * (N + lags) %/% 3] quarterly_sums;
  int<lower=0> measurement_errors_num;
  array[M] int measurement_errors_binary;
  vector<lower=0>[measurement_errors_num] measurement_errors_sigma;
  vector[2] measurement_errors_hyper_prior;
  int<lower=0> measurement_errors_hyper_num;
  int<lower=0> measurement_errors_size;
  int<lower=0, upper=1> parallel_likelihood;
  int<lower=0, upper=1> minnesota_parallel;
  vector[M] data_transformation_mean;
  vector[M] data_transformation_scale;
  int<lower=0, upper=1> hyperbolic_transformation;
  array[M] int<lower=0, upper=1> no_cross_shrinkage_vec;
}
parameters {
  vector<lower=-1, upper=1>[sgt_len_vec[1]] lambda;
  vector<lower=log(p_q_mins[1] + 0.001)>[sgt_len_vec[2]] log_p;
  vector<lower=log(p_q_mins[2] + 0.001)>[sgt_len_vec[3]] log_q;
  vector[M*M-B_zero_res_num-B_pos_res_num-B_neg_res_num] B_par;
  vector<lower=0>[B_pos_res_num] B_pos;
  vector<upper=0>[B_neg_res_num] B_neg;
  vector[M*include_constant] constant;
  vector[M*M*lags-A_zero_res_num] A_par;
  vector<lower=hyper_shrinkage_lim[1], upper=hyper_shrinkage_lim[2]>[hyper_len_vec[1]] hyper_shrinkage;
  vector<lower=0, upper=10>[hyper_len_vec[2]] hyper_lags;
  vector<lower=0, upper=1>[hyper_len_vec[3]] hyper_ownlags;
  vector[soc_free] hyper_soc;
  vector[dio_free] hyper_dio;
  vector<lower=0>[include_hyper_B] hyper_B;
  vector<lower=0>[measurement_errors_hyper_num] hyper_measurement_errors;
  array[(include_garch_groups == 1 ? garch_group_num : M*include_garch)] simplex[(garch_dependence == 1 ? M + 2 : 3)] garch_param;
  array[include_vol_breaks *  M] simplex[vol_breaks_num + 1] relative_vol_raw;
  vector[missing_data_count] missing_data_raw;
  vector[measurement_errors_size] measurement_errors;
  array[number_of_quarterly * (N + lags) %/% 3] simplex[3] monthly_raw;
}
transformed parameters {
  vector[M*include_minnesota] reduced_scales;
  matrix[M, M] B;
  matrix[M*lags, M] A;
  matrix[(include_garch_groups == 1 ? garch_group_num : M*include_garch), (include_garch_groups == 1 ? garch_group_num : M*include_garch)] garch_D;
  vector[(include_garch_groups == 1 ? garch_group_num : M*include_garch)] garch_c;
  vector[(include_garch_groups == 1 ? garch_group_num : M*include_garch)] garch_C;
  matrix<lower=0, upper=(vol_breaks_num+1)>[include_vol_breaks * M, vol_breaks_num + 1] relative_vol;
  vector[missing_data_count] missing_data;
  if(hyperbolic_transformation == 1) {
    missing_data = (exp(missing_data_raw) - exp(-missing_data_raw)) / 2;
  } else {
    missing_data = missing_data_raw;
  }
  for(i in 1:M) {
    for(j in 1:M) {
      if (B_zero_res[(j-1)*M + i, 1] == 1) {
        B[i, j] = 0;
      } else if (B_pos_res[(j-1)*M + i, 1] == 1) {
        B[i, j] = B_pos[B_pos_res[(j-1)*M + i, 2]];
      } else if (B_neg_res[(j-1)*M + i, 1] == 1) {
        B[i, j] = B_neg[B_neg_res[(j-1)*M + i, 2]];
      } else {
        B[i, j] = B_par[(j-1)*M + i - B_zero_res[(j-1)*M + i, 2] - B_pos_res[(j-1)*M + i, 2] - B_neg_res[(j-1)*M + i, 2]];
      }
    }
  }
  {
    int count = 0;
    for(j in 1:M) {
      for(i in 1:(M*lags)) {
        if (A_zero_res[(j-1)*M*lags + i, 1] == 1) {
          A[i, j] = 0;
        } else {
          count += 1;
          A[i, j] = A_par[count];
        }
      }
    }
  }
  if (include_minnesota == 1) {
    if (B_inverse == 0) {
      if (fix_moments == 2) {
        for(i in 1:M) reduced_scales[i] = sqrt(sum(row(B, i) .* row(B, i)));
      } else {
        for(i in 1:M) reduced_scales[i] = sum(fabs(row(B, i)));
      }
    } else {
      for(i in 1:M) reduced_scales[i] = 1.0;
    }
  }
  if (include_garch == 1) {
    if (garch_dependence == 1) {
      for(i in 1:M) {
        garch_c[i] = garch_param[i][1];
        garch_C[i] = garch_param[i][2];
        for(j in 1:M) {
          garch_D[i,j] = garch_param[i][j + 2];
        }
      }
    } else {
      if (include_garch_groups == 0) {
        for(i in 1:M) {
          garch_c[i] = garch_param[i][1];
          garch_C[i] = garch_param[i][2];
          for(j in 1:M) {
            if(i == j) {
              garch_D[i,j] = garch_param[i][3];
            } else {
              garch_D[i,j] = 0.0;
            }
          }
        }
      } else {
        for(i in 1:garch_group_num) {
          garch_c[i] = garch_param[i][1];
          garch_C[i] = garch_param[i][2];
          for(j in 1:garch_group_num) {
            if(i == j) {
              garch_D[i,j] = garch_param[i][3];
            } else {
              garch_D[i,j] = 0.0;
            }
          }
        }
      }
    }
  }
  if (include_vol_breaks == 1) {
    for(i in 1:M) {
      for(j in 1:(vol_breaks_num + 1)) {
        relative_vol[i,j] = relative_vol_raw[i][j] * (vol_breaks_num + 1);
      }
    }
  }
}
model {
  matrix[N, M] xA;
  matrix[N, M] shocks;
  matrix[(N == 0 ? N : N - 1), M] shocks_short;
  int grainsize = 1;
  array[M] int col_indices;
  matrix[M, 3] gamma;
  matrix[N, M] cmat;
  vector[M] cvec;
  array[3] real hypers;
  vector[M*M] Bvec;
  vector[(include_garch_groups == 1 ? garch_group_num : M * include_garch)] garch_v;
  vector[include_vol_breaks * M] vol_v;
  int vol_break_count = 1;
  int parallel_vol_eval;
  matrix[N + lags, M + M*lags] yx_raw;
  matrix[N + lags, M] y_raw_new;
  matrix[N, M] yy;
  matrix[N, M*lags] xx;

  if (N > 0) {

    ///////////////////////////////////////////////////
    // Autoregressive part and missing data business //
    ///////////////////////////////////////////////////

    if (include_constant == 1) {
      cvec = constant;
    } else {
      cvec = rep_vector(0, M);
    }
    cmat = rep_vector(1, N) * cvec';
    if(missing_data_count == 0 && number_of_quarterly == 0 && measurement_errors_num == 0) {
      if (lags > 0) {
        xA = x * A;
      } else {
        xA = rep_matrix(0, N, M);
      }
      if (B_inverse == 0) shocks = (y - cmat - xA) / B';
      if (B_inverse == 1) shocks = (y - cmat - xA) * B';
    } else {
      if(number_of_quarterly > 0) {
        y_raw_new = rebuild_y_raw(y_raw, quarterly_binary, quarterly_sums, monthly_raw, data_transformation_mean, data_transformation_scale, M, (N + lags) %/% 3);
      } else {
        y_raw_new = y_raw;
      }
      if(measurement_errors_num > 0) {
        y_raw_new = subtract_measurement_errors(y_raw_new, measurement_errors_binary, measurement_errors, M, N + lags);
      }
      yx_raw = build_yx(y_raw_new, lags, missing_data, missing_data_location, missing_data_count, N, M);
      yy = yx_raw[(lags + 1):(N + lags), 1:M];
      xx = yx_raw[(lags + 1):(N + lags), (M + 1):(M + M*lags)];
      if (lags > 0) {
        xA = xx * A;
      } else {
        xA = rep_matrix(0, N, M);
      }
      if (B_inverse == 0) shocks = (yy - cmat - xA) / B';
      if (B_inverse == 1) shocks = (yy - cmat - xA) * B';
    }

    ////////////////////////////////////
    // Conditional heteroskedasticity //
    ////////////////////////////////////

    if (include_garch > 0) {

      // Initial values for the shock volatility process (GARCH)
      if (include_garch_groups == 0) garch_v = rep_vector(1, M); else garch_v = rep_vector(1, garch_group_num);

      // Dependent shock volatility processes (efficient parallelization not possible)
      if (garch_dependence == 1) {
        parallel_vol_eval = 0;
        for(i in 2:N) {
          if(garch_eta_form == 0) garch_v = sqrt(garch_c + garch_C .* (garch_v .* garch_v) + garch_D * ( (shocks[i-1,:] .* garch_v') .* (shocks[i-1,:] .* garch_v') )' );
          if(garch_eta_form == 1) garch_v = sqrt(garch_c + garch_C .* (garch_v .* garch_v) + garch_D * (shocks[i-1,:] .* shocks[i-1,:])' );
          shocks[i,:] = shocks[i,:] ./ garch_v';
          target += -sum(log(fabs(garch_v)));
        }

      } else {

        // Independent shock volatility processes (efficient parallelization possible)
        if(include_garch_groups == 0) {
          if(parallel_likelihood == 1) {
            parallel_vol_eval = 1;
          } else {
            parallel_vol_eval = 0;
            for(j in 1:M) {
              for(i in 2:N) {
                if(garch_eta_form == 0) garch_v[j] = sqrt(garch_param[j][1] + garch_param[j][2] * garch_v[j]^2 + garch_param[j][3] * (shocks[i-1,j] * garch_v[j])^2);
                if(garch_eta_form == 1) garch_v[j] = sqrt(garch_param[j][1] + garch_param[j][2] * garch_v[j]^2 + garch_param[j][3] * shocks[i-1,j]^2);
                shocks[i,j] = shocks[i,j] / garch_v[j];
                target += -log(fabs(garch_v[j]));
              }
            }
          }

        // Grouped shock volatility processes (efficient parallelization not implemented)
        } else {
          parallel_vol_eval = 0;
          for(i in 2:N) {
            for(g in 1:garch_group_num) {
              int garch_group_size = 0;
              for(j in 1:M) if(garch_group_mat[j, g] == 1) garch_group_size += 1;
              vector[garch_group_size] last_shocks;
              int count = 1;
              for(j in 1:M) {
                if(garch_group_mat[j, g] == 1) {
                  if(garch_eta_form == 0) last_shocks[count] = shocks[(i-1), j] * garch_v[g];
                  if(garch_eta_form == 1) last_shocks[count] = shocks[(i-1), j];
                  count += 1;
                }
              }
              garch_v[g] = sqrt(garch_param[g][1] + garch_param[g][2] * garch_v[g]^2 + garch_param[g][3] * mean(last_shocks .* last_shocks));
              count = 1;
              for(j in 1:M) {
                if(garch_group_mat[j, g] == 1) {
                  shocks[i,j] = shocks[i,j] / garch_v[g];
                  target += -log(fabs(garch_v[g]));
                  count += 1;
                }
              }
            }
          }
        }
      }
    }

    if(include_vol_breaks == 1) {
      for(i in 1:N) {
        if(vol_breaks[i] == 1) vol_break_count += 1;
        vol_v = relative_vol[:,vol_break_count];
        shocks[i,:] = shocks[i,:] ./ vol_v';
        target += -sum(log(fabs(vol_v)));
      }
    }

    //////////////////////////
    // Log determinant of B //
    //////////////////////////

    if (include_garch == 1) {
      if (B_inverse == 0) target += -(N - 1) * log_determinant(B);
      if (B_inverse == 1) target += (N - 1) * log_determinant(B);
    } else {
      if (B_inverse == 0) target += -N * log_determinant(B);
      if (B_inverse == 1) target += N * log_determinant(B);
    }

  } // if (N > 0) ends

  ///////////////////////////////////
  // Shock distribution parameters //
  ///////////////////////////////////

  for(i in 1:M) {
    if(sgt_len_vec[1] == M) {
      gamma[i, 1] = lambda[i];
    } else {
      gamma[i, 1] = sgt_fixed[1];
    }
    if(sgt_len_vec[2] == M) {
      gamma[i, 2] = log_p[i];
    } else {
      gamma[i, 2] = sgt_fixed[2];
    }
    if(sgt_len_vec[3] == M) {
      gamma[i, 3] = log_q[i];
    } else {
      gamma[i, 3] = sgt_fixed[3];
    }
  }

  //////////////////////////////
  // Shock density evaluation // (i.e. the bulk of likelihood evaluation)
  //////////////////////////////

  if(N > 0) {
    if(include_garch == 1) if(N > 0) shocks_short = shocks[2:N, 1:M];
    if (parallel_likelihood == 0) {
      for(i in 1:M) {
        if(include_garch == 1) {
          target += sgt_lpdf(col(shocks_short, i) | N-1, gamma[i, 1], gamma[i, 2], gamma[i, 3], fix_moments);
        } else {
          target += sgt_lpdf(col(shocks, i) | N, gamma[i, 1], gamma[i, 2], gamma[i, 3], fix_moments);
        }
      }
    } else {
      for(i in 1:M) col_indices[i] = i;
      if(parallel_vol_eval == 1) {
        target += reduce_sum(likelihood_vol_partial_sum, col_indices, grainsize, N, shocks, gamma, fix_moments,
                             garch_param, garch_v, garch_eta_form);
      } else {
        target += reduce_sum(likelihood_partial_sum, col_indices, grainsize, N, shocks, gamma, fix_moments);
      }
    }
  }

  ////////////
  // Priors //
  ////////////

  for(i in 1:M) {
    if (sgt_len_vec[1] == M) target += beta_lpdf((gamma[i, 1] / 2) + 0.5 | lambda_prior[1], lambda_prior[2]);
    if (p_q_prior_shift[1] == 0 && p_q_prior_shift[2] == 0) {
      if (sgt_len_vec[2] == M && p_prior[2] != positive_infinity()) target += normal_lpdf(gamma[i, 2] | p_prior[1], p_prior[2]);
      if (sgt_len_vec[3] == M && q_prior[2] != positive_infinity()) target += normal_lpdf(gamma[i, 3] | q_prior[1], q_prior[2]);
    } else {
      if (sgt_len_vec[2] == M && p_prior[2] != positive_infinity()) target += lognormal_lpdf(exp(gamma[i, 2]) - p_q_prior_shift[1] | p_prior[1], p_prior[2]) + gamma[i, 2];
      if (sgt_len_vec[3] == M && q_prior[2] != positive_infinity()) target += lognormal_lpdf(exp(gamma[i, 3]) - p_q_prior_shift[2] | q_prior[1], q_prior[2]) + gamma[i, 3];
    }
  }
  if (include_minnesota == 1) {
    hypers = hyper_fixed;
    if (hyper_len_vec[1] == 1) hypers[1] = hyper_shrinkage[1];
    if (hyper_len_vec[2] == 1) hypers[2] = hyper_lags[1];
    if (hyper_len_vec[3] == 1) hypers[3] = hyper_ownlags[1];
    if (parallel_likelihood == 0 || minnesota_parallel == 0) {
      for(h in 1:lags) {
        for(i in 1:M) {
          for(j in 1:M) {
            if(A_zero_res[(j-1)*M*lags + (h-1)*M + i, 1] != 1) {
              if (i == j) {
                if (h == 1) {
                  target += normal_lpdf(A[(h-1)*M + i, j] | minnesota_means[i], hypers[1]);
                } else {
                  target += normal_lpdf(A[(h-1)*M + i, j] | 0, hypers[1] / (h^hypers[2]));
                }
              } else {
                if(no_cross_shrinkage_vec[i] + no_cross_shrinkage_vec[j] == 2) {
                  target += normal_lpdf(A[(h-1)*M + i, j] | 0, (hypers[1] * reduced_scales[j]) / ((h^hypers[2]) * reduced_scales[i]) );
                } else {
                  target += normal_lpdf(A[(h-1)*M + i, j] | 0, (hypers[1] * hypers[3] * reduced_scales[j]) / ((h^hypers[2]) * reduced_scales[i]) );
                }
              }
            }
          }
        }
      }
    } else {
      target += reduce_sum(minnesota_partial_sum, col_indices, grainsize, lags, M, A, minnesota_means, hypers, reduced_scales, no_cross_shrinkage_vec, A_zero_res);
    }
    if (include_minnesota == 1 && hyper_len_vec[1] == 1 && hyper_shrinkage_prior[2] != positive_infinity()) target += lognormal_lpdf(hypers[1] | hyper_shrinkage_prior[1], hyper_shrinkage_prior[2]);
    if (include_minnesota == 1 && hyper_len_vec[2] == 1 && hyper_lags_prior[2] != positive_infinity()) target += lognormal_lpdf(hypers[2] | hyper_lags_prior[1], hyper_lags_prior[2]);
    if (include_minnesota == 1 && hyper_len_vec[3] == 1) target += beta_lpdf(hypers[3] | hyper_ownlags_prior[1], hyper_ownlags_prior[2]);
    if (include_soc == 1) {
      real soc_val;
      if (soc_dio_dep[1] == 1) soc_val = soc * hypers[1]; else soc_val = soc;
      if(soc_free == 1) {
        if (soc_dio_dep[1] == 1) soc_val = exp(hyper_soc[1]) * hypers[1]; else soc_val = exp(hyper_soc[1]);
        if (soc_prior[2] != positive_infinity()) target += normal_lpdf(hyper_soc | soc_prior[1], soc_prior[2]);
      }
      for(j in 1:M) {
        real sum_of_coefs = 0;
        for(h in 1:lags) {
          sum_of_coefs += A[(h-1)*M + j, j];
        }
        target += normal_lpdf(sum_of_coefs | minnesota_means[j], soc_val);
      }
    }
    if (include_dio == 1) {
      real dio_val;
      if (soc_dio_dep[2] == 1) dio_val = dio * hypers[1]; else dio_val = dio;
      if(dio_free == 1) {
        if (soc_dio_dep[2] == 1) dio_val = exp(hyper_dio[1]) * hypers[1]; else dio_val = exp(hyper_dio[1]);
        if (dio_prior[2] != positive_infinity()) target += normal_lpdf(hyper_dio | dio_prior[1], dio_prior[2]);
      }
      for(j in 1:M) {
        real sum_of_all_coefs = sum(col(A, j));
        if (include_constant == 1) sum_of_all_coefs += constant[j];
        target += normal_lpdf(sum_of_all_coefs | minnesota_means[j], dio_val);
      }
    }
  }
  if (include_B_prior == 1) {
    for(i in 1:M) Bvec[((i-1)*M+1):(i*M)] = col(B, i);
    if (B_prior_cov_diagonal == 1) {
      int prior_elast_count = 1;
      for(i in 1:(M*M)) {
        if (B_prior_cov[i,i] != positive_infinity() || prior_elast_binary[i] == 1) {
          if (prior_elast_binary[i] == 1) {
            real denumerator = Bvec[prior_elast_ind[prior_elast_count, 2]];
            real elasticity = Bvec[i] / denumerator;
            target += student_t_lpdf( elasticity | prior_elast_par[prior_elast_count, 3], prior_elast_par[prior_elast_count, 1], prior_elast_par[prior_elast_count, 2]);
            target += -log(fabs(denumerator));
            prior_elast_count += 1;
          } else {
            if (include_hyper_B == 1 && B_prior_cov[i,i] == -99) {
              if (B_zero_res[i, 1] == 0) target += normal_lpdf(Bvec[i] | B_prior_mean[i], hyper_B[1] );
            } else {
              if (B_pos_res[i, 1] == 1) {
                target += lognormal_lpdf(Bvec[i] | B_prior_cov[i,i] + log(B_prior_mean[i]), sqrt(B_prior_cov[i,i])); // Log-normal with mode equal to 'B_prior_mean'
              } else if (B_neg_res[i, 1] == 1) {
                target += lognormal_lpdf(-Bvec[i] | B_prior_cov[i,i] + log(-B_prior_mean[i]), sqrt(B_prior_cov[i,i])); // Log-normal with mode equal to 'B_prior_mean'
              } else {
                if (B_zero_res[i, 1] == 0) target += normal_lpdf(Bvec[i] | B_prior_mean[i], (sqrt(B_prior_cov[i,i])) );
              }
            }
          }
        }
      }
    } else {
      // ... (to be implemented?)
    }
  }
  if (include_hyper_B == 1) target += lognormal_lpdf(hyper_B | hyper_B_prior[1], hyper_B_prior[2]);
  if (constant_prior[2] != positive_infinity() && include_constant == 1) {
    for(i in 1:M) target += normal_lpdf(constant[i] | constant_prior[1], constant_prior[2]);
  }
  if (include_garch == 1) {
    if(include_garch_groups == 0) {
      for(i in 1:M) target += dirichlet_lpdf(garch_param[i] | garch_prior);
    } else {
      for(i in 1:garch_group_num) target += dirichlet_lpdf(garch_param[i] | garch_prior);
    }
  }
  if (include_vol_breaks == 1) {
    for(i in 1:M) target += dirichlet_lpdf(relative_vol_raw[i] | vol_breaks_prior);
  }
  if (hyperbolic_transformation == 1 && missing_data_count > 0) {
    for(i in 1:missing_data_count) {
      target += (exp(missing_data_raw[i]) + exp(-missing_data_raw[i])) / 2;
    }
  }
  if (measurement_errors_num > 0 && N > 0) {
    int me_count = 0;
    int me_sigma_count = 0;
    int me_hyper_count = 0;
    for(j in 1:M) {
      if(measurement_errors_binary[j] == 1) {
        me_sigma_count += 1;
        if(measurement_errors_sigma[me_sigma_count] == 0) {
          me_hyper_count += 1;
          target += lognormal_lpdf(hyper_measurement_errors[me_hyper_count] | measurement_errors_hyper_prior[1], measurement_errors_hyper_prior[2]);
        }
        for(i in 1:(N + lags)) {
          if(y_raw[i,j] != positive_infinity()) {
            me_count += 1;
            if(measurement_errors_sigma[me_sigma_count] == 0) {
              target += normal_lpdf(measurement_errors[me_count] | 0, hyper_measurement_errors[me_hyper_count]);
            } else {
              target += normal_lpdf(measurement_errors[me_count] | 0, measurement_errors_sigma[me_sigma_count]);
            }
          }
        }
      }
    }
  }
}
generated quantities {
  vector<lower=p_q_mins[1]>[sgt_len_vec[2]] p;
  vector<lower=p_q_mins[2]>[sgt_len_vec[3]] q;
  vector<lower=0>[soc_free] hyper_soc_exp;
  vector<lower=0>[dio_free] hyper_dio_exp;
  vector<lower=fix_moments>[pq_len] pq;
  vector[missing_data_count] missing_data_real;
  matrix[M, M] B_inv;
  matrix[N + lags, number_of_quarterly] yq;
  vector[prior_elast_num] elasticity;
  matrix[N, M] shocks;
  matrix[N, M] residuals;
  matrix[N, (include_garch_groups == 1 ? garch_group_num : M)] volatility;
  if(include_garch_groups == 1) {
    volatility = rep_matrix(1.0, N, garch_group_num);
  } else {
    volatility = rep_matrix(1.0, N, M);
  }
  if (soc_free == 1) hyper_soc_exp = exp(hyper_soc);
  if (dio_free == 1) hyper_dio_exp = exp(hyper_dio);
  if (sgt_len_vec[2] == M) {
    p = exp(log_p);
  }
  if (sgt_len_vec[3] == M) {
    q = exp(log_q);
  }
  if (sgt_len_vec[2] == M && sgt_len_vec[3] == M) pq = p .* q;
  if (sgt_len_vec[2] == 0 && sgt_len_vec[3] == M) pq = exp(sgt_fixed[2]) * q;
  if (sgt_len_vec[2] == M && sgt_len_vec[3] == 0) pq = p * exp(sgt_fixed[3]);
  if (sgt_len_vec[2] == 0 && sgt_len_vec[3] == 0) pq = rep_vector(exp(sgt_fixed[2]) * exp(sgt_fixed[3]), pq_len);
  B_inv = inverse(B);
  if (missing_data_count > 0) {
    matrix[N + lags, M] y_new = y_raw;
    for(i in 1:missing_data_count) {
      y_new[missing_data_location[i, 1], missing_data_location[i, 2]] = y_new[(missing_data_location[i, 1] - 1), missing_data_location[i, 2]] + missing_data[i];
    }
    for(i in 1:missing_data_count) {
      int y_col = missing_data_location[i, 2];
      int y_row = missing_data_location[i, 1];
      missing_data_real[i] = (y_new[y_row, y_col] * data_transformation_scale[y_col]) + data_transformation_mean[y_col];
    }
  }
  if (number_of_quarterly > 0) {
    matrix[N + lags, M] y_raw_new;
    int count = 1;
    y_raw_new = rebuild_y_raw(y_raw, quarterly_binary, quarterly_sums, monthly_raw, data_transformation_mean, data_transformation_scale, M, (N + lags) %/% 3);
    for(j in 1:M) {
      if(quarterly_binary[j] == 1) {
        yq[:, count] = y_raw_new[:, j];
        count = count + 1;
      }
    }
  }
  if (prior_elast_num > 0) {
    vector[M*M] Bvec;
    for(i in 1:M) Bvec[((i-1)*M+1):(i*M)] = col(B, i);
    int prior_elast_count = 1;
    for(i in 1:(M*M)) {
      if (prior_elast_binary[i] == 1) {
        real denumerator = Bvec[prior_elast_ind[prior_elast_count, 2]];
        elasticity[prior_elast_count] = Bvec[i] / denumerator;
        prior_elast_count += 1;
      }
    }
  }
  if (N > 0) {
    matrix[N, M] xA;
    matrix[N, M] cmat;
    vector[M] cvec;
    vector[(include_garch_groups == 1 ? garch_group_num : M * include_garch)] garch_v;
    vector[include_vol_breaks * M] vol_v;
    int vol_break_count = 1;
    matrix[N + lags, M + M*lags] yx_raw;
    matrix[N + lags, M] y_raw_new;
    matrix[N, M] yy;
    matrix[N, M*lags] xx;
    if (include_constant == 1) {
      cvec = constant;
    } else {
      cvec = rep_vector(0, M);
    }
    cmat = rep_vector(1, N) * cvec';
    if(missing_data_count == 0 && number_of_quarterly == 0) {
      if (lags > 0) {
        xA = x * A;
      } else {
        xA = rep_matrix(0, N, M);
      }
      residuals = (y - cmat - xA);
      if (B_inverse == 0) shocks = (y - cmat - xA) / B';
      if (B_inverse == 1) shocks = (y - cmat - xA) * B';
    } else {
      if(number_of_quarterly > 0) {
        y_raw_new = rebuild_y_raw(y_raw, quarterly_binary, quarterly_sums, monthly_raw, data_transformation_mean, data_transformation_scale, M, (N + lags) %/% 3);
      } else {
        y_raw_new = y_raw;
      }
      if(measurement_errors_num > 0) {
        y_raw_new = subtract_measurement_errors(y_raw_new, measurement_errors_binary, measurement_errors, M, N + lags);
      }
      yx_raw = build_yx(y_raw_new, lags, missing_data, missing_data_location, missing_data_count, N, M);
      yy = yx_raw[(lags + 1):(N + lags), 1:M];
      xx = yx_raw[(lags + 1):(N + lags), (M + 1):(M + M*lags)];
      if (lags > 0) {
        xA = xx * A;
      } else {
        xA = rep_matrix(0, N, M);
      }
      residuals = (yy - cmat - xA);
      if (B_inverse == 0) shocks = (yy - cmat - xA) / B';
      if (B_inverse == 1) shocks = (yy - cmat - xA) * B';
    }
    if (include_garch > 0) {
      if (include_garch_groups == 0) garch_v = rep_vector(1, M); else garch_v = rep_vector(1, garch_group_num);
      volatility[1,:] = garch_v';
      if (garch_dependence == 1) {
        for(i in 2:N) {
          if(garch_eta_form == 0) garch_v = sqrt(garch_c + garch_C .* (garch_v .* garch_v) + garch_D * ( (shocks[i-1,:] .* garch_v') .* (shocks[i-1,:] .* garch_v') )' );
          if(garch_eta_form == 1) garch_v = sqrt(garch_c + garch_C .* (garch_v .* garch_v) + garch_D * (shocks[i-1,:] .* shocks[i-1,:])' );
          shocks[i,:] = shocks[i,:] ./ garch_v';
          volatility[i,:] = garch_v';
        }
      } else {
        if(include_garch_groups == 0) {
          for(j in 1:M) {
              for(i in 2:N) {
                if(garch_eta_form == 0) garch_v[j] = sqrt(garch_param[j][1] + garch_param[j][2] * garch_v[j]^2 + garch_param[j][3] * (shocks[i-1,j] * garch_v[j])^2);
                if(garch_eta_form == 1) garch_v[j] = sqrt(garch_param[j][1] + garch_param[j][2] * garch_v[j]^2 + garch_param[j][3] * shocks[i-1,j]^2);
                shocks[i,j] = shocks[i,j] / garch_v[j];
                volatility[i,j] = garch_v[j];
              }
            }
        } else {
          for(i in 2:N) {
            for(g in 1:garch_group_num) {
              int garch_group_size = 0;
              for(j in 1:M) if(garch_group_mat[j, g] == 1) garch_group_size += 1;
              vector[garch_group_size] last_shocks;
              int count = 1;
              for(j in 1:M) {
                if(garch_group_mat[j, g] == 1) {
                  if(garch_eta_form == 0) last_shocks[count] = shocks[(i-1), j] * garch_v[g];
                  if(garch_eta_form == 1) last_shocks[count] = shocks[(i-1), j];
                  count += 1;
                }
              }
              garch_v[g] = sqrt(garch_param[g][1] + garch_param[g][2] * garch_v[g]^2 + garch_param[g][3] * mean(last_shocks .* last_shocks));
              volatility[i,g] = garch_v[g];
              count = 1;
              for(j in 1:M) {
                if(garch_group_mat[j, g] == 1) {
                  shocks[i,j] = shocks[i,j] / garch_v[g];
                  count += 1;
                }
              }
            }
          }
        }
      }
    }
    if(include_vol_breaks == 1) {
      for(i in 1:N) {
        if(vol_breaks[i] == 1) vol_break_count += 1;
        vol_v = relative_vol[:,vol_break_count];
        shocks[i,:] = shocks[i,:] ./ vol_v';
        volatility[i,:] = vol_v';
      }
    }
  }
}



