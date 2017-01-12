function setKinetics(w, left, right)
% SETKINETICS  Set the surface reaction mechanisms on a wall.
% setKinetics(w, left, right)
% Use class ReactorSurface instead. To be removed after Cantera 2.3.
%
% :param w:
%     Instance of class :mat:func:`Wall`
% :param left:
%     Instance of class :mat:func:`Kinetics` (or another object
%     derived from Kinetics) to be used as the
%     surface kinetics for the left side of the wall. Typically
%     an instance of class :mat:func:`Interface`
% :param right:
%     Instance of class :mat:func:`Kinetics` (or another object
%     derived from Kinetics) to be used as the
%     surface kinetics for the right side of the wall. Typically
%     an instance of class :mat:func:`Interface`
%

warning('This function is deprecated, and will be removed after Cantera 2.3. Use class ReactorSurface instead.');

ileft = 0;
iright = 0;
if isa(left, 'Kinetics')
    ileft = kinetics_hndl(left);
end

if isa(right,'Kinetics')
    iright = kinetics_hndl(right);
end
wallmethods(12, wall_hndl(w), ileft, iright);
