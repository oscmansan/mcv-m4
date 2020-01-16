function X_trian = triangulate(x1, x2, P1, P2, imsize)

    %Normalization step:
    %Scaling and translation so that both pixel coord. are in interval [-1,1]

    nx = imsize(1);
    ny = imsize(2);

    H = [2/nx,  0,   -1;
          0,  2/ny,  -1;
          0,    0,    1];

    x1_h = [x1;1];
    x2_h = [x2;1]; 

    x1 = H*x1_h;
    x2 = H*x2_h;
    P1 = H*P1;
    P2 = H*P2;


    %Convert to homogeneous coordinates
    x1 = [x1(1)/x1(3), x1(2)/x1(3), 1];
    x2 = [x2(1)/x2(3), x2(2)/x2(3), 1];

    %Build the constraint matrix
    A = [x1(1)*P1(3,:) - P1(1,:);
         x1(2)*P1(3,:) - P1(2,:);
         x2(1)*P2(3,:) - P2(1,:);
         x2(2)*P2(3,:) - P2(2,:)];

    [~,~,V] =svd(A,0); % use the economy decomposition

    % Take the last column of the transposed of V, that's the singular
    % vector with the lowest singular value.
    X_trian = V(:,4);

     
