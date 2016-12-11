function [u] = G7_sol_Alternative_Laplace_Equation_Axb(f, dom2Inp, param)
%this code is not intended to be efficient. 

[ni, nj]=size(f);

%We add the ghost boundaries (for the boundary conditions)
f_ext = zeros(ni+4, nj+4);
f_ext(3:end-2, 3:end-2) = f;
dom2Inp_ext = zeros(ni+4, nj+4);
dom2Inp_ext (3:end-2, 3:end-2) = dom2Inp;

%Store memory for the A matrix and the b vector    
nPixels =(ni+2)*(nj+2); %Number of pixels

%We will create A sparse, this is the number of nonzero positions

%idx_Ai: Vector for the nonZero i index of matrix A
%idx_Aj: Vector for the nonZero j index of matrix A
%a_ij: Vector for the value at position ij of matrix A
values = 5*sum(sum(dom2Inp)) + 2*(nj + ni + 4) + ni*nj - sum(sum(dom2Inp)); 

idx_Ai = zeros(values, 1);
idx_Aj = zeros(values, 1);
a_ij = zeros(values, 1);


b = zeros(nPixels,1);

%Vector counter
idx = 1;

%North side boundary conditions
i = 1;
for j = 1:nj+4
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+4) + i;
    
    
    %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
    %vector b
    idx_Ai(idx) = p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx = idx + 1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p + 2;
    a_ij(idx) = -1;   
    idx=idx+1;
            
    b(p) = 0;
end

i = 2;
for j = 1:nj+4
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+4) + i;
    
    
    %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
    %vector b
    idx_Ai(idx) = p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx = idx + 1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p + 1;
    a_ij(idx) = -1;   
    idx=idx+1;
            
    b(p) = 0;
end

%South side boundary conditions
i = ni+3;
for j = 1:nj+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;
    
    %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
    %vector b
    idx_Ai(idx) = p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx = idx+1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p-1;
    a_ij(idx) = -1;   
    idx = idx + 1;
            
    b(p) = 0;
    
end

i = ni+4;
for j = 1:nj+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;
    
    %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
    %vector b
    idx_Ai(idx) = p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx = idx+1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p-2;
    a_ij(idx) = -1;   
    idx = idx + 1;
            
    b(p) = 0;
    
end

%West side boundary conditions
j=1;
for i=3:ni+2
% for i=1:ni+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;

    %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
    %vector b
    idx_Ai(idx) = p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx = idx+1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p+ni+4;
    a_ij(idx) = -1;   
    idx = idx + 1;
            
    b(p) = 0;
    
    
end

j = 2;
for i=3:ni+2
% for i=1:ni+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;

    %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
    %vector b
    idx_Ai(idx) = p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx = idx+1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p+ni+3;
    a_ij(idx) = -1;   
    idx = idx + 1;
            
    b(p) = 0;
    
    
end

%East side boundary conditions
j=nj+3;
for i=1:ni+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;
    
    %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
    %vector b
    idx_Ai(idx) = p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx = idx+1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p-(ni+4);
    a_ij(idx) = -1;   
    idx = idx + 1;
            
    b(p) = 0;
    
end

j=nj+4;
for i=1:ni+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;
    
    %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
    %vector b
    idx_Ai(idx) = p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx = idx+1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p-(ni+3);
    a_ij(idx) = -1;   
    idx = idx + 1;
            
    b(p) = 0;
    
end

%Inner points
for j=3:nj+2
    for i=3:ni+2
     
        %from image matrix (i,j) coordinates to vectorial (p) coordinate
        p = (j-1)*(ni+4)+i;
                                            
        if (dom2Inp_ext(i,j)==1) %If we have to inpaint this pixel
            
            %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
            %vector b
            
            idx_Ai(idx) = p; 
            idx_Aj(idx) = p-(2*(ni+4)); 
            a_ij(idx) = 1/(2*(param.hj)^2);
            idx = idx+1;

            idx_Ai(idx) = p; 
            idx_Aj(idx) = p-(ni+4); 
            a_ij(idx) = -1/(2*(param.hj)^2);
            idx = idx+1;
    
            idx_Ai(idx) = p; 
            idx_Aj(idx) = p-2; 
            a_ij(idx) = 1/(2*(param.hi)^2);
            idx = idx+1;
 
            idx_Ai(idx) = p; 
            idx_Aj(idx) = p-1; 
            a_ij(idx) = -1/(2*(param.hi)^2);
            idx = idx+1;

            idx_Ai(idx) = p; 
            idx_Aj(idx) = p; 
            a_ij(idx) = (-1/2)*((1/(param.hi)^2)+(1/(param.hj)^2));
            idx = idx+1;
            
            idx_Ai(idx) = p; 
            idx_Aj(idx) = p+1; 
            a_ij(idx) = 1/(2*(param.hi)^2);
            idx = idx+1;
            
            idx_Ai(idx) = p; 
            idx_Aj(idx) = p+ni+4; 
            a_ij(idx) = 1/(2*(param.hj)^2);
            idx = idx+1;

            b(p) = 0;
            
    
        else %we do not have to inpaint this pixel 
            
            %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
            %vector b
            idx_Ai(idx) = p; 
            idx_Aj(idx) = p; 
            a_ij(idx) = 1;
            idx = idx+1;

            b(p) = f_ext(i,j);
            
        end       
    end
end
    %A is a sparse matrix, so for memory requirements we create a sparse
    %matrix
    %TO COMPLETE 7
    A=sparse(idx_Ai, idx_Aj, a_ij, (ni+2)*(nj+2), (ni+2)*(nj+2)); %??? and ???? is the size of matrix A
    
    %Solve the sistem of equations
    x=mldivide(A,b);
    
    %From vector to matrix
    u_ext= reshape(x, ni+2, nj+2);
    
    %Eliminate the ghost boundaries
    u=full(u_ext(2:end-1, 2:end-1));
    
