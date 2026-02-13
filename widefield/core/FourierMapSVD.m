function  [ComplexMap, AbsMap, AngleMap] = FourierMapSVD(U, V, TimeVec, ff)

Tensor = svdFrameReconstruct(U, V);

[ComplexMap, AbsMap, AngleMap] = FourierMap(Tensor, TimeVec, ff);