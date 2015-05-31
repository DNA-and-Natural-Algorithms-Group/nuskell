#
# The section 6, figure 8 "cooperative hybridization" CRN implementation from Thachuk, Winfree, Soloveichik, LNCS xxx, 2015
# 
# Coded by Erik Winfree (winfree@caltech.edu)
#
# Clamps are not implemented.
#
# In order to enumerate cooperative hybridization without considering all toehold release steps to be slow,
# we use two toehold sizes, and just the "right-side" toeholds are "slow".  Our solution involves enumeration
# using the options          --release-cutoff-1-1 4 --release-cutoff-1-n 8 --k-fast 0.1
# and thus the verifying compiler invocation will be
#
# ./verify ts/soloveichik_cooperative.ts crn/crn4z.crn --release-cutoff-1-1 4 --release-cutoff-1-n 8 --k-fast 0.1

class formal(s) = "? ? ? d1a d1b d1c d2a d2b d2c"
                | ". . .  .   .   .   .   .   . "
    where {
        d1a = unique(5); 
        d1b = unique(5);
        d1c = unique(7);
        d2a = unique(5); 
        d2b = unique(5);
        d2c = unique(7) };

module bifuel1(A,B,X,Y,Q) = 
       "A1b A1c A2a A2b A2c B1a B1b B1c B2a B2b Q1a Q1b Q1c + B2c* B2b* B2a* B1c* B1b* B1a* A2c* A2b* A2a* A1c* A1b* A1a*"
     | " (   (   (   (   (   (   (   (   (   (   .   .   .  +  .    )    )    )    )    )    )    )    )    )    )    .  "
    where {
    	A1a = A.d1a; A1b = A.d1b; A1c = A.d1c;
    	A2a = A.d2a; A2b = A.d2b; A2c = A.d2c;
	B1a = B.d1a; B1b = B.d1b; B1c = B.d1c;
	B2a = B.d2a; B2b = B.d2b; B2c = B.d2c;
	Q1a = Q.d1a; Q1b = Q.d1b; Q1c = Q.d1c };

module bifuel2(A,B,X,Y,Q) = 
       "A2b A2c B1a B1b B1c X1a X1b X1c + B2a B2b Q1a Q1b Q1c Y1a Y1b Y1c + Q1c* Q1b* Q1a* B2b* B2a* B1c* B1b* B1a* A2c* A2b* A2a*"
     | " (   (   (   (   (   .   .   .  +  (   (   (   (   (   .   .   .  +  )    )    )    )    )    )    )    )    )    )    .  "
    where {
    	A1a = A.d1a; A1b = A.d1b; A1c = A.d1c;
    	A2a = A.d2a; A2b = A.d2b; A2c = A.d2c;
	B1a = B.d1a; B1b = B.d1b; B1c = B.d1c;
	B2a = B.d2a; B2b = B.d2b; B2c = B.d2c;
	X1a = X.d1a; X1b = X.d1b; X1c = X.d1c;
	Y1a = Y.d1a; Y1b = Y.d1b; Y1c = Y.d1c;
	Q1a = Q.d1a; Q1b = Q.d1b; Q1c = Q.d1c };

module bifuel3(A,B,X,Y,Q) = 
       "B1b B1c X1a X1b X1c X2a X2b X2c + X1c* X1b* X1a* B1c* B1b* B1a*"
     | " (   (   (   (   (   .   .   .  +  )    )    )    )    )    .  "
    where {
	B1a = B.d1a; B1b = B.d1b; B1c = B.d1c;
	X1a = X.d1a; X1b = X.d1b; X1c = X.d1c;
	X2a = X.d2a; X2b = X.d2b; X2c = X.d2c };

module bifuel4(A,B,X,Y,Q) = 
       "Q1b Q1c Y1a Y1b Y1c Y2a Y2b Y2c + Y1c* Y1b* Y1a* Q1c* Q1b* Q1a*"
     | " (   (   (   (   (   .   .   .  +  )    )    )    )    )    .  "
    where {
	Y1a = Y.d1a; Y1b = Y.d1b; Y1c = Y.d1c;
	Y2a = Y.d2a; Y2b = Y.d2b; Y2c = Y.d2c;
	Q1a = Q.d1a; Q1b = Q.d1b; Q1c = Q.d1c };

# This molecule is never actually created; it is just used to pass the Q domain conveniently.
macro commonQ(r) = "d1a d1b d1c" | ". . ."
    where {
        d1a = unique(5);
        d1b = unique(5);
        d1c = unique(7) };

module bimolecular(r) = infty(bifuel1(reactantA,reactantB,productX,productY,Q)) 
       		      + infty(bifuel2(reactantA,reactantB,productX,productY,Q)) 
		      + infty(bifuel3(reactantA,reactantB,productX,productY,Q)) 
       		      + infty(bifuel4(reactantA,reactantB,productX,productY,Q)) 
    where {
    	Q = commonQ(r);
        reactantA = r.reactants[0];
	reactantB = r.reactants[1];
	productX = r.products[0];
	productY = r.products[1] };

module rxn(r) =
    if len(r.reactants) == 2 and len(r.products) == 2 then
        bimolecular(r)
    else r[0];
# this last line is a hack to cause a crash when the target CRN is invalid for this scheme

module main(crn) = sum(map(rxn, crn))
    where
        crn = irrev_reactions(crn)
