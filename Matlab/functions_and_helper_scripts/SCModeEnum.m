classdef SCModeEnum < Simulink.IntEnumType
  enumeration
    STABILIZE(0)
    SLEW(1)
    DETUMBLE(2)
    DESATURATE(3)
  end
end