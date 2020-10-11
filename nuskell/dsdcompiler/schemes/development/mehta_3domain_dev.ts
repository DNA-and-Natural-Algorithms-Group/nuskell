#
# A "three domain" scheme by Jenish Mehta, 2015.  
# Each signal species is a single long domain flanked by a toehold on each side.
# 
# Coded by Erik Winfree (winfree@caltech.edu)
#
# This should be modified to take care of all different arities.  Right now it is just 2-to-2.
# So it can only be run on crn4, crn17, and reaction1... reaction14
#
# There appears to be an undesirably 4-way branch migration pathway (not yet investigated).
# Also, the waste from the irreversible pathway is not inert, which causes problems for
# pathway decomposition, and therefore only bisimulation deems it to be correct.
#
# ./verify ts/mehta_3domain.ts crn/crn4.crn --bisimulation --ignore-branch-4way



global toe_s = short();
global toe_t = short();
global toe_w = short();
global irrev = long();

class formal(s) = "toe_s branch toe_t"
                | "  .     .      .  "
    where {
        branch = long() }; 

module revbigate_forward(sA,sB,sC,sD) = 
       "A s + B w + s C t + s D t + t* D* t* C* w* B* s* A* s*"
     | "( ( + ( ( + . ( ( + . ( ( + )  )  )  )  )  )  )  )  . "
    where {
        s = toe_s; t = toe_t; w = toe_w;
    	A = sA.branch; B = sB.branch; C = sC.branch; D = sD.branch };

module revbigate_backward(sA,sB,sC,sD) = 
       "s A t + s B t + w C + t D + t* D* t* C* w* B* s* A* s*"
     | "( ( . + ( ( . + ( ( + ( ( + .  )  )  )  )  )  )  )  ) "
    where {
        s = toe_s; t = toe_t; w = toe_w;
    	A = sA.branch; B = sB.branch; C = sC.branch; D = sD.branch };

module revbihelpers(sA,sB,sC,sD) = 
       [ "A s" | ". .", "B w" | ". .", "w C" | ". .", "t D" | ". ." ]
    where {
        s = toe_s; t = toe_t; w = toe_w;
    	A = sA.branch; B = sB.branch; C = sC.branch; D = sD.branch };

module revbimolecular(r) = 
       infty(revbigate_forward(reactantA,reactantB,productC,productD)) 
     + infty(revbigate_backward(reactantA,reactantB,productC,productD)) 
     + sum(map(infty, revbihelpers(reactantA,reactantB,productC,productD)))
    where {
        reactantA = r.reactants[0];
	reactantB = r.reactants[1];
	productC = r.products[0];
	productD = r.products[1] };

module irrevbigate_forward(sA,sB,sC,sD) = 
       "A s + B w + I + I w + s C t + s D t + t* D* t* C* w* I* I* w* B* s* A* s*"
     | "( ( + ( ( + ( + ( ( + . ( ( + . ( ( + )  )  )  )  )  )  )  )  )  )  )  . "
    where {
        s = toe_s; t = toe_t; w = toe_w; I = irrev;
    	A = sA.branch; B = sB.branch; C = sC.branch; D = sD.branch };

module irrevbihelpers(sA,sB,sC,sD) = 
       [ "A s" | ". .", "B w" | ". .", "w C" | ". .", "t D" | ". .", "w I I" | ". . ." ]
    where {
        s = toe_s; t = toe_t; w = toe_w; I = irrev;
    	A = sA.branch; B = sB.branch; C = sC.branch; D = sD.branch };

module irrevbimolecular(r) = 
       infty(irrevbigate_forward(reactantA,reactantB,productC,productD)) 
     + sum(map(infty, irrevbihelpers(reactantA,reactantB,productC,productD)))
    where {
        reactantA = r.reactants[0];
	reactantB = r.reactants[1];
	productC = r.products[0];
	productD = r.products[1] };

module rxn(r) =
    if len(r.reactants) == 2 and len(r.products) == 2 and r.reversible then
        revbimolecular(r)
    elseif len(r.reactants) == 2 and len(r.products) == 2 then
        irrevbimolecular(r)
    else r[0];
# this last line is a hack to cause a crash when the target CRN is invalid for this scheme.
# actually, this scheme can be easily generalized to higher and lower arity reactions. and should be...

module main(crn) = sum(map(rxn, crn))
    where
        crn = rev_reactions(crn)


# Here it is in Visual DSD, as per March 11 email to Jenish Mehta

# def Fuel = 20
# 
# def signal(N,A) =
#   ( N * <s^ A t^> )
# 
# (* reversible reaction *)
# def birev(A,B,C,D) =
#   ( constant Fuel * {s^*}[A s^]:[B w^]:<s^>[C t^]:<s^>[D t^]
#   | constant Fuel * <A s^>
#   | constant Fuel * <B w^>
#   | constant Fuel * <w^ C>
#   | constant Fuel * <t^ D>
#   | constant Fuel * [s^ A]<t^>:[s^ B]<t^>:[w^ C]:[t^ D]{t^*}
#   )
# 
# (* irreversible reaction *)
# def biirrev(A,B,C,D) =
#   ( constant Fuel * {s^*}[A s^]:[B w^]:[i]:[i w^]:<s^>[C t^]:<s^>[D t^]
#   | constant Fuel * <A s^>
#   | constant Fuel * <B w^>
#   | constant Fuel * <w^ C>
#   | constant Fuel * <t^ D>
#   | constant Fuel * <w^ i i>
#   )
# 
# 
# ( birev(A,B,C,D)
# | signal(5,A)
# | signal(10,B)
# | biirrev(P,Q,R,S)
# | signal(5,P)
# | signal(10,Q)
# )
# 