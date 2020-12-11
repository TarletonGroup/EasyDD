function Q = findRotationMatrixAtoB(a, b)
    % Find the rotation matrix that brings vector a into b.
    % https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
    a = a / norm(a);
    b = b / norm(b);
    v = cross(a, b);
    c = dot(a, b);
    V = [0 -v(3) v(2);
        v(3) 0 -v(1);
        -v(2) v(1) 0];
    Q = eye(3) + V + (V^2) / (1 + c);
end
