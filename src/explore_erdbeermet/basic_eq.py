from sympy import symbols
a1,a2=symbols("a1 a2")
vx0,ux0,xy0,vy0,uy0, uv0=symbols("vx0 ux0 xy0 vy0 uy0 uv0")
dx1, dy1, dz1= symbols("dx dy dz")
# x,y:z
vx1=vx0+dx1
ux1=ux0+dx1
ux1=ux0+dx1
uy1=uy0+dy1
vy1=vy0+dy1
uz1=a1*ux0+(1-a1)*uy0+dz1
vz1=a1*vx0+(1-a1)*vy0+dz1
uz1=a1*ux0+(1-a1)*uy0+dz1
xy1=xy0+dx1+dy1
xz1=(1-a1)*xy0+dz1+dx1
yz1=a1*xy0+dy1+dz1
def calc_a(uz,vy,vz,uy,ux,vx):
    return(((uz+vy)-(vz+uy))/((ux+vy)-(vx+uy)))
# x,y:q
xq2,yq2,uq2,vq2,yq2,zq2,dq2,dx2,dy2,uv1=symbols("xq2 yq2 uq2 vq2 yq2 zq2 dq2 dx2 dy2 uv1")
xq2=(1-a2)*xy1+dq2+dx2
yq=a2*xy1+dy2+dq2
uq2=a2*ux1+(1-a2)*uy1+dq2
vq2=a2*vx1+(1-a2)*vy1+dq2
zq2=a2*xz1+(1-a2)*yz1+dq2

xy2=xy1+dx2+dy2
ux2=ux1+dx2
vx2=vx1+dx2
xz2=xz1+dx2
uy2=uy1+dy2
vy2=vy1+dy2
yz2=yz1+dy2
uz2=uz1
xz2=xz1+dx2
yz2=yz1+dy2
vz2=vz1
uv1=uv0
