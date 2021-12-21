
n_stages = 9
a = [[0]*(n_stages+1) for _ in range(n_stages + 1)]
c = [0] * (n_stages+1)
b = [0]* (n_stages+1)
bh = [0]*(n_stages+1)

c[1] =  0
c[2] =  3/50
c[3] =  1439/15000
c[4] =  1439/10000
c[5] =  4973/10000
c[6] =  389/400
c[7] =  1999/2000
c[8] =  1
c[9] =  1
# 
# ********************************************************
#  COUPLING COEFFICIENTS
#  ---------------------
#for c[1] =  0
#
# for c[2] =  3/50
a[2][1] =  3/50
#  
# for c[3] =  1439/15000
a[3][1] =  519479/27000000
a[3][2] =  2070721/27000000
#  
# for c[4] =  1439/10000
a[4][1] =  1439/40000
a[4][2] =  0
a[4][3] =  4317/40000
#  
# for c[5] =  4973/10000
a[5][1] =  109225017611/82828840000
a[5][2] =  0
a[5][3] = -417627820623/82828840000
a[5][4] =  43699198143/10353605000
#  
# for c[6] =  389/400
a[6][1] = -8036815292643907349452552172369/191934985946683241245914401600
a[6][2] =  0
a[6][3] =  246134619571490020064824665/1543816496655405117602368
a[6][4] = -13880495956885686234074067279/113663489566254201783474344
a[6][5] =  755005057777788994734129/136485922925633667082436
# 
# for c[7] =  1999/2000
a[7][1] = -1663299841566102097180506666498880934230261/30558424506156170307020957791311384232000
a[7][2] =  0
a[7][3] =  130838124195285491799043628811093033/631862949514135618861563657970240
a[7][4] = -3287100453856023634160618787153901962873/20724314915376755629135711026851409200
a[7][5] =  2771826790140332140865242520369241/396438716042723436917079980147600
a[7][6] = -1799166916139193/96743806114007800
# 
# for c[8] =  1
a[8][1] = -832144750039369683895428386437986853923637763/15222974550069600748763651844667619945204887
a[8][2] =  0
a[8][3] =  818622075710363565982285196611368750/3936576237903728151856072395343129
a[8][4] = -9818985165491658464841194581385463434793741875/61642597962658994069869370923196463581866011
a[8][5] =  31796692141848558720425711042548134769375/4530254033500045975557858016006308628092
a[8][6] = -14064542118843830075/766928748264306853644
a[8][7] = -1424670304836288125/2782839104764768088217
# 
# for c[9] =  1
a[9][1] =  382735282417/11129397249634
a[9][2] =  0
a[9][3] =  0
a[9][4] =  5535620703125000/21434089949505429
a[9][5] =  13867056347656250/32943296570459319
a[9][6] =  626271188750/142160006043
a[9][7] = -51160788125000/289890548217
a[9][8] =  163193540017/946795234
# 
# ********************************************************
#  High order weights c[9] =  1
#  (This is also the propagating stage 9, with  a[9,i] =  b[i].)
#  --------------------------------------------------------
# 
b[1] =  382735282417/11129397249634
b[2] =  0
b[3] =  0
b[4] =  5535620703125000/21434089949505429
b[5] =  13867056347656250/32943296570459319
b[6] =  626271188750/142160006043
b[7] = -51160788125000/289890548217
b[8] =  163193540017/946795234
b[9] =  0
# 
# ********************************************************
#  Low order weights with  c[extra] =  1
#  --------------------------------------------------
# 
bh[1] =  124310637869885675646798613/2890072468789466426596827670
bh[2] =  0
bh[3] =  0
bh[4] =  265863151737164990361330921875/1113197271463372303940319369579
bh[5] =  3075493557174030806536302953125/6843749922042323876546949699876
bh[6] =  67798000008733879813263055/29532792147666737550036372
bh[7] = -1099436585155390846238326375/15055706496446408859196167
bh[8] =  26171252653086373181571802/368794478890732346033505
bh[9] =  1/30

A = [[0]* n_stages for _ in range(n_stages)]
C = [0] * n_stages
B = [0] * n_stages
B_HAT = [0] * n_stages

print(a)

for i in range(len(A)):
    for j in range(len(A[i])):
        A[i][j] = a[i + 1][j + 1]

for i in range(len(C)):
    C[i] = c[i + 1]

for i in range(len(B)):
    B[i] = b[i + 1]

for i in range(len(B_HAT)):
    B_HAT[i] = bh[i + 1]

print("Checking if the values were set up properly")
for i in range(len(A)):
    print("A", A[i])
    print("a", a[i+1])
print("=====================================================")

print(B)
print(b)
print("===================================================")

print(B_HAT)
print(bh)
print("===================================================")

print(C)
print(c)
print("===================================================")


print("getting the pairs")

print("A =", A)
print("C =", C)
print("B =", B)
print("B_HAT =", B_HAT)