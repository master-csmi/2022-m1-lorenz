#include <parareal/utils.hpp>

Vector<double> lorenz(double /*t*/, Vector<double> X, double* gamma){
    Vector<double> sol(X.cols());
    sol << gamma[0] * (X[1]-X[0]), X[0] * (gamma[2]-X[2])-X[1],
         X[0]*X[1]-gamma[1]*X[2];

    return sol;
}

Matrix RK4(Vector<double> X0, double dt, double t0, int nb_t, 
        Vector<double> prob(double, Vector<double>, double*), double* gamma){

    int dim = static_cast<int>(X0.cols());
    Matrix X(1,dim); X << X0;

    double t = t0;
    Vector<double> K1(dim), K2(dim), K3(dim), K4(dim);
    Vector<double> X_prec(dim);
    
    // while ( (t+dt)<=T or std::abs(t+dt-T)<1e-6){
    // while ( (t+dt)<=T ){ 
    for(int i=1; i<nb_t; i++){
        X_prec = X.bottomRows<1>();
        K1=prob(t, X_prec, gamma);
        K2=prob(t+dt/2., X_prec + 1./2. * K1 * dt, gamma);
        K3=prob(t+dt/2., X_prec + 1./2. * K2 * dt, gamma);
        K4=prob(t+dt, X_prec+ K3 * dt, gamma);

        X.conservativeResize(X.rows()+1, X.cols());
        X.row(X.rows()-1) = X_prec + dt/6.* (K1+2.*K2+2.*K3+K4);

        t+=dt;
    }

    return X;
}

Vector<double> compute_times(double t0, double T, double dt_G, double dt_F, int n_proc,
        Vector<int>* tab_nb_t_G_p, Vector<int>* tab_nb_t_F_p){

    int nb_t_G = static_cast<int>((T-t0)/dt_G);
    int nb_t_F = static_cast<int>((T-t0)/dt_F);
    int rapport = static_cast<int>(dt_G/dt_F);

    // to compute the number of coarse time step for each process
    for(int i=0; i<n_proc; i++){
        (*tab_nb_t_G_p)(i) = nb_t_G/n_proc;
    }

    int reste = nb_t_G-tab_nb_t_G_p->sum();
    int index=0;
    for(int i=0; i<reste; i++){
        if(index==n_proc-1)
            index=0;
        (*tab_nb_t_G_p)(index)++;
        index++;
    }
    
    // to compute the number of fine time step for each process
    // with final coarse time == final fine time (on each interval)
    for(int i=0; i<n_proc-1; i++){
        (*tab_nb_t_F_p)(i) = (*tab_nb_t_G_p)(i) * rapport;
    }
    (*tab_nb_t_F_p)(n_proc-1) = nb_t_F-tab_nb_t_F_p->head(n_proc-1).sum(); 

    // to avoid problems with size
    assert(nb_t_G==tab_nb_t_G_p->sum());
    assert(nb_t_F==tab_nb_t_F_p->sum());

    // to compute the t_j
    Vector<double> times(n_proc+1);
    times(0) = t0;
    for(int i=1; i<n_proc+1; i++){
        times(i) = times(i-1) + (*tab_nb_t_G_p)(i-1) * dt_G;
    }

    return times;
}

bool sol_converge(Matrix X0_k, Matrix X0_knext, double eps){
    Matrix diff = X0_knext-X0_k;
    return diff.cwiseAbs().maxCoeff() < eps;
}