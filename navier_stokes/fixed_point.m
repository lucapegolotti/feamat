function u = fixed_point(u_old,dofsu,fespace_u,M,H,b,dt,tol,kmax)

count = 0;
u = u_old;

nodes1 = 1:dofsu;
nodes2 = dofsu+1:2*dofsu;

C = assemble_convective_term(fespace_u,u);
C = apply_dirichlet_bc_matrix(C,fespace_u,0);

Hit = H;
Hit(nodes1,nodes1) = Hit(nodes1,nodes1) + C;
Hit(nodes2,nodes2) = Hit(nodes2,nodes2) + C;

err = norm(M*(u-u_old)/dt + Hit*u - b);

while (err > tol && count < kmax)
    count = count + 1;
    
    u = (M+dt*Hit)\(M*u_old+dt*b);
    
    C = assemble_convective_term(fespace_u,u);
    C = apply_dirichlet_bc_matrix(C,fespace_u,0);

    Hit = H;
    Hit(nodes1,nodes1) = Hit(nodes1,nodes1) + C;
    Hit(nodes2,nodes2) = Hit(nodes2,nodes2) + C;
    
    err = norm(M*(u-u_old)/dt + Hit*u - b);
    
    disp(['Fixed point iteration number ', num2str(count),', error = ', num2str(err)]);
end

if (count >= kmax)
    error('We reached the maximum number of iterations for fixed point!');
end