#include <parareal/utils.hpp>

Vector lorenz(double /*t*/, Vector X, int /*dim*/, double* gamma){
    Vector sol(X.cols());
    sol << gamma[0] * (X[1]-X[0]), X[0] * (gamma[2]-X[2])-X[1],
         X[0]*X[1]-gamma[1]*X[2];

    return sol;
}

Matrix RK4(Vector X0, double dt, double t0, double T, Vector prob(double, Vector, 
        int, double*), double* gamma){

    int dim = X0.cols();
    Matrix X(1,dim); X << X0;

    double t = t0;
    Vector K1(dim), K2(dim), K3(dim), K4(dim);
    Vector X_prec(dim);
    
    while ( (t+dt)<=T or std::abs(t+dt-T)<1e-6){ 
        X_prec = X.bottomRows<1>();
        K1=prob(t, X_prec, dim, gamma);
        K2=prob(t+dt/2., X_prec + 1./2. * K1 * dt, dim, gamma);
        K3=prob(t+dt/2., X_prec + 1./2. * K2 * dt, dim, gamma);
        K4=prob(t+dt, X_prec+ K3 * dt, dim, gamma);

        X.conservativeResize(X.rows()+1, X.cols());
        X.row(X.rows()-1) = X_prec + dt/6.* (K1+2.*K2+2.*K3+K4);

        t+=dt;
    }

    return X;
}

Vector compute_times(double t0, double T, double dt_G, int P){
    // time between t_j and t_{j+1}
    double dt_P = (T-t0)/P;
    // nb points
    int nb_pts = dt_P/dt_G;

    // t_j exact
    Vector times_exact(P+1);
    times_exact[0] = t0;
    for(int j=1; j<P+1; j++){
        times_exact[j] = times_exact[j-1] + dt_P;
    }
    times_exact[P]=T;

    // t_j approach (to be a multiple of dt_G)
    double exact, approach;
    Vector times(P+1);
    times[0] = t0;
    for(int j=1; j<P+1; j++){
        exact=times_exact[j];
        approach=times[j-1]+dt_G*nb_pts;
        if(exact-approach<=abs(exact-(approach+dt_G))) { // lower rounding
            times[j] = times[j-1] + dt_G*nb_pts;
        }
        else { // upper rounding
            times[j] = times[j-1] + dt_G*(nb_pts+1);
        }
    }
    times[P]=T;

    return times;
}

bool sol_converge(Matrix X0_k, Matrix X0_knext, double eps){
    Matrix diff = X0_knext-X0_k;
    return diff.cwiseAbs().maxCoeff() < eps;
}