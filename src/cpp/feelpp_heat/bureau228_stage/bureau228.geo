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
Hide "*";
//+
Show {
  Point{10}; Point{12}; Point{14}; Point{16}; Point{18}; Point{20}; Point{22}; Point{24}; Curve{28}; Curve{30}; Curve{32}; Curve{34}; Curve{36}; Curve{38}; Curve{60}; Curve{62}; Curve{64}; Curve{66}; Surface{160}; Surface{164}; Surface{168}; 
}
//+
Hide "*";
//+
Show {
  Point{10}; Point{12}; Point{14}; Point{16}; Point{18}; Point{20}; Point{22}; Point{24}; Curve{28}; Curve{30}; Curve{32}; Curve{34}; Curve{36}; Curve{38}; Curve{60}; Curve{62}; Curve{64}; Curve{66}; Surface{160}; Surface{164}; Surface{168}; 
}
//+
Hide "*";
//+
Show {
  Point{10}; Point{12}; Point{14}; Point{16}; Point{18}; Point{20}; Point{22}; Point{24}; Curve{28}; Curve{30}; Curve{32}; Curve{34}; Curve{36}; Curve{38}; Curve{60}; Curve{62}; Curve{64}; Curve{66}; Surface{160}; Surface{164}; Surface{168}; 
}
//+
Hide "*";
//+
Show {
  Point{3}; Point{5}; Point{11}; Point{13}; Point{19}; Point{21}; Curve{23}; Curve{29}; Curve{35}; Curve{53}; Curve{55}; Curve{61}; Curve{63}; Surface{150}; Surface{162}; 
}
//+
Hide "*";
//+
Show {
  Point{1}; Point{2}; Point{3}; Point{4}; Point{5}; Point{6}; Point{7}; Point{8}; Curve{1}; Curve{2}; Curve{3}; Curve{4}; Curve{21}; Curve{22}; Curve{23}; Curve{24}; Curve{25}; Curve{26}; Surface{84}; Surface{86}; Surface{88}; 
}