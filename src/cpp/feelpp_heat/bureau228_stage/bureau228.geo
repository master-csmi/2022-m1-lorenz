h=1;

ws={0,2.4,3.4,3.7};
npw=#ws[];
nlw=npw-1;
ls={0,5.0};
npl=#ls[];
nll=npl-1;
zs={0,0.7,1.7,1.9,2.3};
npz=#zs[];
nlz=npz-1;

pts={};
For z In {0: npz-1}
    For w In {0: npw-1}
        For l In {0: npl-1}
            pts[]+={newp};
            Point(newp) = {ls[l],ws[w],zs[z],h};
        EndFor
    EndFor
EndFor

lns={};
For z In {0: npz-1}
    For w In {0: npw-1}
        For l In {0: npl-2}
            lns[]+={newl};
            Line(newl) = {pts[z*npw*npl+w*npl+l],pts[z*npw*npl+w*npl+l+1]};
        EndFor
    EndFor
EndFor
For z In {0: npz-1}
    For w In {0: npw-2}
        For l In {0: npl-1}
            lns[]+={newl};
            Line(newl) = {pts[z*npw*npl+w*npl+l],pts[z*npw*npl+(w+1)*npl+l]};
        EndFor
    EndFor
EndFor
For z In {0: npz-2}
    For w In {0: npw-1}
        For l In {0: npl-1}
            lns[]+={newl};
            Line(newl) = {pts[z*npw*npl+w*npl+l],pts[(z+1)*npw*npl+w*npl+l]};
        EndFor
    EndFor
EndFor

nl1=npz*npw*nll;
nl2=nl1+npz*nlw*npl;
sfs={};
For z In {0: npz-1}
    For w In {0: nlw-1}
        For l In {0: nll-1}
            Line Loop(newll) = {-lns[z*nll*npw+w*nll+l],lns[nl1+z*npl*nlw+w*npl+l],lns[z*nll*npw+w*nll+l+1],-lns[nl1+z*npl*nlw+w*npl+1]};
            sfs[]+={news};
            Surface(news) = {newll-1};
        EndFor
    EndFor
EndFor
For z In {0: nlz-1}
    For w In {0: npw-1}
        For l In {0: nll-1}
            Line Loop(newll) = {-lns[z*nll*npw+w*nll+l],lns[nl2+z*npl*npw+w*npl+l],lns[(z+1)*nll*npw+w*nll+l],-lns[nl2+z*npl*npw+w*npl+l+1]};
            sfs[]+={news};
            Surface(news) = {newll-1};
        EndFor
    EndFor
EndFor
For z In {0: nlz-1}
    For w In {0: nlw-1}
        For l In {0: npl-1}
            Line Loop(newll) = {-lns[nl1+z*nlw*npl+w*npl+l],lns[nl2+z*npl*npw+w*npl+l],lns[nl1+(z+1)*nlw*npl+w*npl+l],-lns[nl2+z*npl*npw+(w+1)*npl+l]};
            sfs[]+={news};
            Surface(news) = {newll-1};
        EndFor
    EndFor
EndFor
// For i In {0:#sfs[]-1}
//     Printf("%g", sfs[i]);
// EndFor

ns1=nll*nlw*npz;
ns2=ns1+nll*npw*nlz;
vms={};
For z In {0: nlz-1}
    For w In {0: nlw-1}
        For l In {0: nll-1}
            Surface Loop(newsl) = {sfs[z*nll*nlw+w*nll+l],sfs[(z+1)*nll*nlw+w*nll+l],sfs[nll*nlw*npz+z*nll*npw+w*nll+l],sfs[nll*nlw*npz+z*nll*npw+(w+1)*nll+l],sfs[ns2+z*npl*nlw+w*npl+l],sfs[ns2+z*npl*nlw+w*npl+l+1]};
            vms[]+={newv};
            Volume(newv) = {newsl-1};
        EndFor
    EndFor
EndFor

Physical Volume("Room") = vms[];
s={};
For i In {0:nll*nlw-1}
    s[]+={sfs[i]};
EndFor
Physical Surface("Floor") = s[];
s={};
For i In {0:nll*nlz-1}
    s[]+={sfs[ns1+i*nll*npw]};
EndFor
Physical Surface("RightWall") = s[];
s={};
For i In {0:nll*nlz-1}
    s[]+={sfs[ns1+(i+1)*nll*npw-1]};
EndFor
Physical Surface("LeftWall") = s[];
s1={sfs[ns2+npl],sfs[ns2+npl+nlw*npl]};
Physical Surface("Door") = s1[];
s={};
For i In {0:nlz*nlw-1}
    s[]+={sfs[ns2+i*npl]};
EndFor
s[]-=s1[];
Physical Surface("FrontWall") = s[];
s={};
For i In {0:nlw-1}
    s[]+={sfs[ns2+i*npl+1]};
EndFor
Physical Surface("LowerBackWall") = s[];
s={};
For i In {0:nlw-1}
    s[]+={sfs[ns2+nlw*npl+i*npl+1]};
EndFor
Physical Surface("LowerWindow") = s[];
s={};
For i In {0:nlw-1}
    s[]+={sfs[ns2+2*nlw*npl+i*npl+1]};
EndFor
Physical Surface("UpperBackWall") = s[];
s={};
For i In {0:nlw-1}
    s[]+={sfs[ns2+3*nlw*npl+i*npl+1]};
EndFor
Physical Surface("UpperWindow") = s[];
s={};
For i In {0:nll*nlw-1}
    s[]+={sfs[nll*nlw*(npz-1)+i]};
EndFor
Physical Surface("Ceiling") = s[];
//+
Point(41) = {-12.0,-4.8,7.0, h};
//+
Point(42) = {-12.7,-6.8,7.0, h};
//+
Point(43) = {-10.0,-3.3,7.0, h};
//+
Point(44) = {-12.0,-6.8,6.8, h};
//+
Point(45) = {-10.5,-6.8,7.0, h};
//+
Point(46) = {-10.0,-5.7,7.0, h};
//+
Point(47) = {-13.5,-3.3,7.0, h};
//+
Point(48) = {-13.0,-2.0,7.0, h};
//+
Point(49) = {-13.5,-5.7,7.0, h};
//+
Point(50) = {-11.5,-2.0,7.0, h};

//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{47}; Point{48}; Point{41}; Point{49}; Point{42}; Point{44}; Point{45}; Point{46}; Point{43}; Point{50}; 
}
//+
Translate {-2, 13.5, -6.3} {
  Point{43}; Point{46}; Point{45}; Point{41}; Point{50}; Point{48}; Point{47}; Point{44}; Point{42}; Point{49}; 
}