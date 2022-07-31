"""This is sub-package is about the Schur-Weyl duality; which is a deep
connetions between permutations of tensor indices, integer partitions
of the tensor rank, and the Tensor's irreducible sub-spaces (under rotations-
in ANY number of dimensions: mind blown).

Seriously: the Schur-Weyl duality is so far from obvious. It is NOT simple.

There are the following modules:

monte.py:
=====
The Symmetric Group: Permutations on n-letters.

young.py:
========
Young Tableaux and Diagrams

pasal.py:
=========
From Pascal's Triangle to Interger Partitions.


Here's the short course:

You get a 9-dimensional rank-2 tensor; what are it's irreducible 'parts'?
We know the answer:

>>>T_S = (T.ij + T.ji) / 2   # The 6-dimensional symmetric part

>>>T_A = (T.ij - T.ji) / 2   # The 3-dimensional antisymmetric part


That was easy. Now note that the tensor rank, 2, has two integer
partitions:

>>>2 == 2

>>>2 == 1 + 1


Each of those can be described by a diagram with 2 boxes:

[ ][ ]

and

[ ]
[ ]


and each of those can be filled with the number (0,1) with rows and columns
strictly increasing as follows:

[0][1]

and

[0]
[1]


Now we look at the Symmetric (permutation) group of order 2:

01 --> 01 = e  # the identity, or neutral element
01 --> 10 = p

and we decide which permutations leave the diamgrams with the same numbers
in the rows:

R_s = (e, p)

R_a = (e,)

and likewise for columns:

C_s = (e,)

C_a = (e, p)

and the we pull a rabbit out of the hat: define the diagrams symmetrizer
as the (normalized) product:

S = sgn(c)C * R / 2!   where we include the sign of the permutations from

the column symmetrizer. Hence we get 2 results:


S_s = (e,) * (e, p) / 2 = (e + p) / 2

S_a - (e, -p) * e / 2 = (e - p) /2

Now we apply those to the tensor indices:


>>>T_S = (T.ij + T.ji) / 2

>>>T_A = (T.ij - T.ji) / 2


Wow-- that *was* c o m p l i c a t e d! But: it is the general algorithm
for a tensor of ANY rank in ANY number of dimensions. That is HUGE.

Each invariant subspace is called a Weyl module, so when you take a rank-N
tensor and say:

>>>for tableau, tensor in T.weyl_modules():
...<suite>


The code takes N = T.rank, and finds all the integer partitions:

>>>for p in pascal.P(N).partitions():

and the for each p, which is a list of integers, it makes a diagram:

>>>d = young.Diagram(*p)

and then it finds all the standard fillings for that:


>>>for tableau in d.standard_tableaux():
...<suite>


It takes that tableau and gets:

R = tableau.Row()
C = tableau.Col()

and then

S = tableau.S()

and then apply that to tensor indices-- AND IT WORKS!


Here is a concrete example:


>>>T = EPSILON & EPSILON

>>>N = T.rank  # 6

>>>p = pascal.P(N)

>>>list(p.partitions())
[[6],
 [5, 1],
 [4, 2],
 [3, 3],
 [4, 1, 1],
 [3, 2, 1],
 [2, 2, 2],
 [3, 1, 1, 1],
 [2, 2, 1, 1],
 [2, 1, 1, 1, 1],
 [1, 1, 1, 1, 1, 1]]


Lets pick [2, 2, 2]:


>>>d = young.Diagram(2, 2, 2)
>>>print d
[ ][ ]
[ ][ ]
[ ][ ]


Now pick a filling (we could have found them all with d.standard_tableaux()):

>>>t = d.fill(0, 3, 1, 4, 2, 5)
>>>print t
[i][l]
[j][m]
[k][n]

What this means is that the tensor is symmetric under interchange
of the (0, 3) and (1, 4) and (2, 5) indices, and antisymmetric in
(0, 1, 2) and (3, 4, 5)--which should be expected considering it is the
symmetric outer product of the rank-3 antisymmetric tensor with itself,
but we still have to do the work:


1st get the Row symmetrizer:


>>>list(t.Row())
[, (25), (14), (14)(25), (03), (03)(25), (03)(14), (03)(14)(25)]


The column symmetrizer pairs each permutation with its sign:

>>>list(t.Col())
[(1, ),
 (-1, (45)),
 (-1, (34)),
 (1, (345)),
 (1, (354)),
 (-1, (35)),
 (-1, (12)),
 (1, (12)(45)),
 (1, (12)(34)),
 (-1, (12)(345)),
 (-1, (12)(354)),
 (1, (12)(35)),
 (-1, (01)),
 (1, (01)(45)),
 (1, (01)(34)),
 (-1, (01)(345)),
 (-1, (01)(354)),
 (1, (01)(35)),
 (1, (012)),
 (-1, (012)(45)),
 (-1, (012)(34)),
 (1, (012)(345)),
 (1, (012)(354)),
 (-1, (012)(35)),
 (1, (021)),
 (-1, (021)(45)),
 (-1, (021)(34)),
 (1, (021)(345)),
 (1, (021)(354)),
 (-1, (021)(35)),
 (-1, (02)),
 (1, (02)(45)),
 (1, (02)(34)),
 (-1, (02)(345)),
 (-1, (02)(354)),
 (1, (02)(35))]

Now we take their product (C*R), invoking the distributive property
(there are a lot of combinatins):


>>>list(t.S())
[(1, ),
 (1, (25)),
 (1, (14)),
 (1, (14)(25)),
 (1, (03)),
 (1, (03)(25)),
 (1, (03)(14)),
 (1, (03)(14)(25)),
 (-1, (45)),
 (-1, (254)),
 (-1, (145)),
 (-1, (1425)),
 (-1, (03)(45)),
 (-1, (03)(254)),
 (-1, (03)(145)),
 (-1, (03)(1425)),
 (-1, (34)),
 (-1, (25)(34)),
 (-1, (143)),
 (-1, (143)(25)),
 (-1, (034)),
 (-1, (034)(25)),
 (-1, (0314)),
 (-1, (0314)(25)),
 (1, (345)),
 (1, (2534)),
 (1, (1453)),
 (1, (14253)),
 (1, (0345)),
 (1, (03425)),
 (1, (03145)),
 (1, (031425)),
 (1, (354)),
 (1, (2543)),
 (1, (1435)),
 (1, (14325)),
 (1, (0354)),
 (1, (03254)),
 (1, (03514)),
 (1, (032514)),
 (-1, (35)),
 (-1, (253)),
 (-1, (14)(35)),
 (-1, (14)(253)),
 (-1, (035)),
 (-1, (0325)),
 (-1, (035)(14)),
 (-1, (0325)(14)),
 (-1, (12)),
 (-1, (152)),
 (-1, (124)),
 (-1, (1524)),
 (-1, (03)(12)),
 (-1, (03)(152)),
 (-1, (03)(124)),
 (-1, (03)(1524)),
 (1, (12)(45)),
 (1, (1542)),
 (1, (1245)),
 (1, (15)(24)),
 (1, (03)(12)(45)),
 (1, (03)(1542)),
 (1, (03)(1245)),
 (1, (03)(15)(24)),
 (1, (12)(34)),
 (1, (152)(34)),
 (1, (1243)),
 (1, (15243)),
 (1, (034)(12)),
 (1, (034)(152)),
 (1, (03124)),
 (1, (031524)),
 (-1, (12)(345)),
 (-1, (15342)),
 (-1, (12453)),
 (-1, (153)(24)),
 (-1, (0345)(12)),
 (-1, (034215)),
 (-1, (031245)),
 (-1, (0315)(24)),
 (-1, (12)(354)),
 (-1, (15432)),
 (-1, (12435)),
 (-1, (15)(243)),
 (-1, (0354)(12)),
 (-1, (032154)),
 (-1, (035124)),
 (-1, (0324)(15)),
 (1, (12)(35)),
 (1, (1532)),
 (1, (124)(35)),
 (1, (15324)),
 (1, (035)(12)),
 (1, (03215)),
 (1, (035)(124)),
 (1, (032415)),
 (-1, (01)),
 (-1, (01)(25)),
 (-1, (041)),
 (-1, (041)(25)),
 (-1, (013)),
 (-1, (013)(25)),
 (-1, (0413)),
 (-1, (0413)(25)),
 (1, (01)(45)),
 (1, (01)(254)),
 (1, (0451)),
 (1, (04251)),
 (1, (013)(45)),
 (1, (013)(254)),
 (1, (04513)),
 (1, (042513)),
 (1, (01)(34)),
 (1, (01)(25)(34)),
 (1, (0431)),
 (1, (0431)(25)),
 (1, (0134)),
 (1, (0134)(25)),
 (1, (04)(13)),
 (1, (04)(13)(25)),
 (-1, (01)(345)),
 (-1, (01)(2534)),
 (-1, (04531)),
 (-1, (042531)),
 (-1, (01345)),
 (-1, (013425)),
 (-1, (045)(13)),
 (-1, (0425)(13)),
 (-1, (01)(354)),
 (-1, (01)(2543)),
 (-1, (04351)),
 (-1, (043251)),
 (-1, (01354)),
 (-1, (013254)),
 (-1, (04)(135)),
 (-1, (04)(1325)),
 (1, (01)(35)),
 (1, (01)(253)),
 (1, (041)(35)),
 (1, (041)(253)),
 (1, (0135)),
 (1, (01325)),
 (1, (04135)),
 (1, (041325)),
 (1, (012)),
 (1, (0152)),
 (1, (0412)),
 (1, (04152)),
 (1, (0123)),
 (1, (01523)),
 (1, (04123)),
 (1, (041523)),
 (-1, (012)(45)),
 (-1, (01542)),
 (-1, (04512)),
 (-1, (042)(15)),
 (-1, (0123)(45)),
 (-1, (015423)),
 (-1, (045123)),
 (-1, (0423)(15)),
 (-1, (012)(34)),
 (-1, (0152)(34)),
 (-1, (04312)),
 (-1, (043152)),
 (-1, (01234)),
 (-1, (015234)),
 (-1, (04)(123)),
 (-1, (04)(1523)),
 (1, (012)(345)),
 (1, (015342)),
 (1, (045312)),
 (1, (042)(153)),
 (1, (012345)),
 (1, (015)(234)),
 (1, (045)(123)),
 (1, (042315)),
 (1, (012)(354)),
 (1, (015432)),
 (1, (043512)),
 (1, (0432)(15)),
 (1, (012354)),
 (1, (0154)(23)),
 (1, (04)(1235)),
 (1, (04)(15)(23)),
 (-1, (012)(35)),
 (-1, (01532)),
 (-1, (0412)(35)),
 (-1, (041532)),
 (-1, (01235)),
 (-1, (015)(23)),
 (-1, (041235)),
 (-1, (0415)(23)),
 (1, (021)),
 (1, (0521)),
 (1, (0241)),
 (1, (05241)),
 (1, (0213)),
 (1, (05213)),
 (1, (02413)),
 (1, (052413)),
 (-1, (021)(45)),
 (-1, (05421)),
 (-1, (02451)),
 (-1, (051)(24)),
 (-1, (0213)(45)),
 (-1, (054213)),
 (-1, (024513)),
 (-1, (0513)(24)),
 (-1, (021)(34)),
 (-1, (0521)(34)),
 (-1, (02431)),
 (-1, (052431)),
 (-1, (02134)),
 (-1, (052134)),
 (-1, (024)(13)),
 (-1, (0524)(13)),
 (1, (021)(345)),
 (1, (053421)),
 (1, (024531)),
 (1, (0531)(24)),
 (1, (021345)),
 (1, (05)(1342)),
 (1, (0245)(13)),
 (1, (05)(13)(24)),
 (1, (021)(354)),
 (1, (054321)),
 (1, (024351)),
 (1, (051)(243)),
 (1, (021354)),
 (1, (054)(132)),
 (1, (024)(135)),
 (1, (051324)),
 (-1, (021)(35)),
 (-1, (05321)),
 (-1, (0241)(35)),
 (-1, (053241)),
 (-1, (02135)),
 (-1, (05)(132)),
 (-1, (024135)),
 (-1, (05)(1324)),
 (-1, (02)),
 (-1, (052)),
 (-1, (02)(14)),
 (-1, (052)(14)),
 (-1, (023)),
 (-1, (0523)),
 (-1, (023)(14)),
 (-1, (0523)(14)),
 (1, (02)(45)),
 (1, (0542)),
 (1, (02)(145)),
 (1, (05142)),
 (1, (023)(45)),
 (1, (05423)),
 (1, (023)(145)),
 (1, (051423)),
 (1, (02)(34)),
 (1, (052)(34)),
 (1, (02)(143)),
 (1, (052)(143)),
 (1, (0234)),
 (1, (05234)),
 (1, (02314)),
 (1, (052314)),
 (-1, (02)(345)),
 (-1, (05342)),
 (-1, (02)(1453)),
 (-1, (053142)),
 (-1, (02345)),
 (-1, (05)(234)),
 (-1, (023145)),
 (-1, (05)(1423)),
 (-1, (02)(354)),
 (-1, (05432)),
 (-1, (02)(1435)),
 (-1, (051432)),
 (-1, (02354)),
 (-1, (054)(23)),
 (-1, (023514)),
 (-1, (0514)(23)),
 (1, (02)(35)),
 (1, (0532)),
 (1, (02)(14)(35)),
 (1, (0532)(14)),
 (1, (0235)),
 (1, (05)(23)),
 (1, (0235)(14)),
 (1, (05)(14)(23))]

To get the symmetrized tensor, one just applies all those permutations
to the tensor indices (via the transpose method), and then sums them
up with the sign included (and a normalization factor that is computed
from the Young diagram via the amazing Hook Length Formula). And a miracle
happens: that symmetrized tensor is closed under rotations!

For curiosity's sake,  there is a method that will convert that
lexigraphically:


>>> t.lexigraphicS()
' (T.ijklmn + T.ijnlmk + T.imkljn + T.imnljk + T.ljkimn + T.ljnimk +
T.lmkijn + T.lmnijk - T.ijklnm - T.ijnlkm - T.imklnj - T.imnlkj - T.ljkinm -
T.ljnikm - T.lmkinj - T.lmnikj - T.ijkmln - T.ijnmlk - T.imkjln - T.imnjlk -
T.ljkmin - T.ljnmik - T.lmkjin - T.lmnjik + T.ijkmnl + T.ijnmkl + T.imkjnl +
T.imnjkl + T.ljkmni + T.ljnmki + T.lmkjni + T.lmnjki + T.ijknlm + T.ijnklm +
T.imknlj + T.imnklj + T.ljknim + T.ljnkim + T.lmknij + T.lmnkij - T.ijknml -
T.ijnkml - T.imknjl - T.imnkjl - T.ljknmi - T.ljnkmi - T.lmknji - T.lmnkji -
T.ikjlmn - T.injlmk - T.ikmljn - T.inmljk - T.lkjimn - T.lnjimk - T.lkmijn -
T.lnmijk + T.ikjlnm + T.injlkm + T.ikmlnj + T.inmlkj + T.lkjinm + T.lnjikm +
T.lkminj + T.lnmikj + T.ikjmln + T.injmlk + T.ikmjln + T.inmjlk + T.lkjmin +
T.lnjmik + T.lkmjin + T.lnmjik - T.ikjmnl - T.injmkl - T.ikmjnl - T.inmjkl -
T.lkjmni - T.lnjmki - T.lkmjni - T.lnmjki - T.ikjnlm - T.injklm - T.ikmnlj -
T.inmklj - T.lkjnim - T.lnjkim - T.lkmnij - T.lnmkij + T.ikjnml + T.injkml +
T.ikmnjl + T.inmkjl + T.lkjnmi + T.lnjkmi + T.lkmnji + T.lnmkji - T.jiklmn -
T.jinlmk - T.mikljn - T.minljk - T.jlkimn - T.jlnimk - T.mlkijn - T.mlnijk +
T.jiklnm + T.jinlkm + T.miklnj + T.minlkj + T.jlkinm + T.jlnikm + T.mlkinj +
T.mlnikj + T.jikmln + T.jinmlk + T.mikjln + T.minjlk + T.jlkmin + T.jlnmik +
T.mlkjin + T.mlnjik - T.jikmnl - T.jinmkl - T.mikjnl - T.minjkl - T.jlkmni -
T.jlnmki - T.mlkjni - T.mlnjki - T.jiknlm - T.jinklm - T.miknlj - T.minklj -
T.jlknim - T.jlnkim - T.mlknij - T.mlnkij + T.jiknml + T.jinkml + T.miknjl +
T.minkjl + T.jlknmi + T.jlnkmi + T.mlknji + T.mlnkji + T.jkilmn + T.jnilmk +
T.mkiljn + T.mniljk + T.jklimn + T.jnlimk + T.mklijn + T.mnlijk - T.jkilnm -
T.jnilkm - T.mkilnj - T.mnilkj - T.jklinm - T.jnlikm - T.mklinj - T.mnlikj -
T.jkimln - T.jnimlk - T.mkijln - T.mnijlk - T.jklmin - T.jnlmik - T.mkljin -
T.mnljik + T.jkimnl + T.jnimkl + T.mkijnl + T.mnijkl + T.jklmni + T.jnlmki +
T.mkljni + T.mnljki + T.jkinlm + T.jniklm + T.mkinlj + T.mniklj + T.jklnim +
T.jnlkim + T.mklnij + T.mnlkij - T.jkinml - T.jnikml - T.mkinjl - T.mnikjl -
T.jklnmi - T.jnlkmi - T.mklnji - T.mnlkji + T.kijlmn + T.nijlmk + T.kimljn +
T.nimljk + T.kljimn + T.nljimk + T.klmijn + T.nlmijk - T.kijlnm - T.nijlkm -
T.kimlnj - T.nimlkj - T.kljinm - T.nljikm - T.klminj - T.nlmikj - T.kijmln -
T.nijmlk - T.kimjln - T.nimjlk - T.kljmin - T.nljmik - T.klmjin - T.nlmjik +
T.kijmnl + T.nijmkl + T.kimjnl + T.nimjkl + T.kljmni + T.nljmki + T.klmjni +
T.nlmjki + T.kijnlm + T.nijklm + T.kimnlj + T.nimklj + T.kljnim + T.nljkim +
T.klmnij + T.nlmkij - T.kijnml - T.nijkml - T.kimnjl - T.nimkjl - T.kljnmi -
T.nljkmi - T.klmnji - T.nlmkji - T.kjilmn - T.njilmk - T.kmiljn - T.nmiljk -
T.kjlimn - T.njlimk - T.kmlijn - T.nmlijk + T.kjilnm + T.njilkm + T.kmilnj +
T.nmilkj + T.kjlinm + T.njlikm + T.kmlinj + T.nmlikj + T.kjimln + T.njimlk +
T.kmijln + T.nmijlk + T.kjlmin + T.njlmik + T.kmljin + T.nmljik - T.kjimnl -
T.njimkl - T.kmijnl - T.nmijkl - T.kjlmni - T.njlmki - T.kmljni - T.nmljki -
T.kjinlm - T.njiklm - T.kminlj - T.nmiklj - T.kjlnim - T.njlkim - T.kmlnij -
T.nmlkij + T.kjinml + T.njikml + T.kminjl + T.nmikjl + T.kjlnmi + T.njlkmi +
T.kmlnji + T.nmlkji) / 720

which can be sent to the eval function, or computed as follows:


>>>print T.symmetric(2, 2, 2)
=============================
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
-----------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
-----------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
=============================
==============================
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
------------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
------------------------------
[0, 0, 0] [0, 0, -1] [0, 1, 0]
[0, 0, 1] [0, 0, 0] [-1, 0, 0]
[0, -1, 0] [1, 0, 0] [0, 0, 0]
==============================
==============================
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
------------------------------
[0, 0, 0] [0, 0, 1] [0, -1, 0]
[0, 0, -1] [0, 0, 0] [1, 0, 0]
[0, 1, 0] [-1, 0, 0] [0, 0, 0]
------------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
==============================
==============================
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
------------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
------------------------------
[0, 0, 0] [0, 0, 1] [0, -1, 0]
[0, 0, -1] [0, 0, 0] [1, 0, 0]
[0, 1, 0] [-1, 0, 0] [0, 0, 0]
==============================
=============================
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
-----------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
-----------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
=============================
==============================
[0, 0, 0] [0, 0, -1] [0, 1, 0]
[0, 0, 1] [0, 0, 0] [-1, 0, 0]
[0, -1, 0] [1, 0, 0] [0, 0, 0]
------------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
------------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
==============================
==============================
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
------------------------------
[0, 0, 0] [0, 0, -1] [0, 1, 0]
[0, 0, 1] [0, 0, 0] [-1, 0, 0]
[0, -1, 0] [1, 0, 0] [0, 0, 0]
------------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
==============================
==============================
[0, 0, 0] [0, 0, 1] [0, -1, 0]
[0, 0, -1] [0, 0, 0] [1, 0, 0]
[0, 1, 0] [-1, 0, 0] [0, 0, 0]
------------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
------------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
==============================
=============================
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
-----------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
-----------------------------
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
[0, 0, 0] [0, 0, 0] [0, 0, 0]
=============================


There a total of 76 isotropies of a rank 6 tensor:

>>>sum([x.dim() for x in itertools.starmap(young.Diagram, p.partitions())])
76

Have fun, and stay safe.
"""
## \namespace geo.metric.schur_weyl The
# <a href="http://en.wikipedia.org/wiki/Schur-Weyl.duality">Schur-Weyl
# Duality, </a> or \f$  S.k \times {\bf GL}(N, \bf{C}) \rightarrow
# (\bf{C}^N)^{\otimes k} \f$

from . import monte
from . import young
from . import pascal
from . import de_bruijn


## This is the crown jewel of the subpackage: yield all symmetries of
# a rank-N tensor (in 1-line of code)- each one is closed under rotations,
# see geo.metric.euclid.euclid.Tensor_.weyl_modules().
# \param tensor_ rank \f$ N \f$ made by euclid.euclid.ranked
# \returns \f$ {\bf T}_{\pi_\lambda} \forall \pi_\lambda \in \Sigma_N \f$
def symmetrize(tensor_):
    """Yield ALL symmetries of a tensor:

    for partish in parition.P(T.rank):
    ....for tab in young.Diagram(*partish).standard_taableaux():
    ........yield T.symmetrize(tab)  # or tab.symmetrize(T)


    that is:
    The rank N can be paritioned P(N) ways (P is the partition function).
    Each parition makes a unique Young diagram (via starmap())
    Each diagram makes 1 or more Standard Young Tableau (whence, chain())
    Each Tableau makes a unique symmetry (Weyl module).



    >>>for p in pascal.P(t.rank).partitions():
    >>>    for tab in young.Diagram(*p).standard_tableaux():
    >>>        yield t.symmetrize(tab)
    """
    from itertools import chain, starmap, imap
    return imap(
        tensor_.symmetrize,
        chain(*(d.standard_tableaux() for d in
                starmap(young.Diagram, pascal.P(tensor_.rank).partitions()))))
