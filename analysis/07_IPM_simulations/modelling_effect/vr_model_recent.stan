data{
    int N;
    int N1;
    int N2;
    int N3;
    int yr_3[N3];
    vector[N2] lagged_2;
    int yr_2[N2];
    vector[N1] lagged_1;
    int yr_1[N1];
    int surv[N1];
    vector[N1] recent_1;
    int z_1[N1];
    vector[N2] z1;
    vector[N2] recent_2;
    int fp[N2];
    int z_2[N2];
    vector[N3] n_seeds;
    int z_3[N3];
    int n_seeds_t[N];
    int n_recr[N];
    vector[N] z_new;
}
parameters{
    real s_int;
    real s_slope;
    real s_recent;
    real g_int;
    real g_slope;
    real g_recent;
    real<lower=0> g_sd;
    real fp_int;
    real fp_slope;
    real seed_int;
    real seed_slope;
    real<lower=0> seed_sd;
    real p_germ;
    real fd_int;
    real<lower=0> fd_sd;
}
model{
    vector[N1] p_s;
    vector[N2] g_mean;
    vector[N2] p_fp;
    vector[N3] seed_mu;
    real seed_p;
    fd_sd ~ exponential( 1 );
    fd_int ~ normal( 0 , 1 );
    z_new ~ normal( fd_int , fd_sd );
    p_germ ~ normal( 0 , 1 );
    seed_p = p_germ;
    seed_p = inv_logit(seed_p);
    n_recr ~ binomial( n_seeds_t , seed_p );
    seed_sd ~ exponential( 1 );
    seed_slope ~ normal( 0 , 1 );
    seed_int ~ normal( 0 , 1 );
    for ( i in 1:N3 ) {
        seed_mu[i] = seed_int + seed_slope * log(z_3[i]);
        seed_mu[i] = exp(seed_mu[i]);
    }
    n_seeds ~ normal( seed_mu , seed_sd );
    fp_slope ~ normal( 0 , 1 );
    fp_int ~ normal( 0 , 1 );
    for ( i in 1:N2 ) {
        p_fp[i] = fp_int + fp_slope * log(z_2[i]);
        p_fp[i] = inv_logit(p_fp[i]);
    }
    fp ~ binomial( 1 , p_fp );
    g_sd ~ exponential( 1 );
    g_recent ~ normal( 0 , 1 );
    g_slope ~ normal( 0 , 1 );
    g_int ~ normal( 0 , 1 );
    for ( i in 1:N2 ) {
        g_mean[i] = g_int + g_slope * log(z_2[i]) + g_recent * recent_2[i];
        g_mean[i] = exp(g_mean[i]);
    }
    z1 ~ normal( g_mean , g_sd );
    s_recent ~ normal( 0 , 1 );
    s_slope ~ normal( 0 , 1 );
    s_int ~ normal( 0 , 1 );
    for ( i in 1:N1 ) {
        p_s[i] = s_int + s_slope * log(z_1[i]) + s_recent * recent_1[i];
        p_s[i] = inv_logit(p_s[i]);
    }
    surv ~ binomial( 1 , p_s );
}
generated quantities{
    vector[N1] log_lik;
    vector[N1] p_s;
    vector[N2] g_mean;
    vector[N2] p_fp;
    vector[N3] seed_mu;
    real seed_p;
    seed_p = p_germ;
    seed_p = inv_logit(seed_p);
    for ( i in 1:N3 ) {
        seed_mu[i] = seed_int + seed_slope * log(z_3[i]);
        seed_mu[i] = exp(seed_mu[i]);
    }
    for ( i in 1:N2 ) {
        p_fp[i] = fp_int + fp_slope * log(z_2[i]);
        p_fp[i] = inv_logit(p_fp[i]);
    }
    for ( i in 1:N2 ) {
        g_mean[i] = g_int + g_slope * log(z_2[i]) + g_recent * recent_2[i];
        g_mean[i] = exp(g_mean[i]);
    }
    for ( i in 1:N1 ) {
        p_s[i] = s_int + s_slope * log(z_1[i]) + s_recent * recent_1[i];
        p_s[i] = inv_logit(p_s[i]);
    }
    for ( i in 1:N1 ) log_lik[i] = binomial_lpmf( surv[i] | 1 , p_s[i] );
}