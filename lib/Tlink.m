function A = Tlink(teta,l,d,alpha)
    % Transformation associated to a link
    % teta -> joint angle (rad)
    % l -> joint length
    % d -> joint displacement
    % alpha -> link torsion angle (rad)
    
    A = rotz(teta)*trans(l,0,d)*rotx(alpha);
end

