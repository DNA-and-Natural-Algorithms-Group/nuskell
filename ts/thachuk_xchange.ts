#
# The section 6, figure 9 "xchange gates" from Thachuk, Winfree, Soloveichik, LNCS xxx, 2015
# 
# Coded by Erik Winfree (winfree@caltech.edu)
#
# Note that there are two "wastes" produced by the standard reaction pathway that aren't formally waste, because they are not inert.
# It is not necessary to start with them present, but they will cause bisimulation and pathway decomposition verification to fail.
#
# Subtle point:  Unlike in figure 9, the unique-to-reaction domain "th" is implemented here as distinct for the "left side" fuel and the "right side" fuel.
# This should also function correctly. In thachuk_xchange_WF.ts, the "th" domain is the same on both sides, as in figure 9.
#
# To avoid combinatorial explosion, this must be invoked along the lines of:
# ./verify ts/thachuk_xchange.ts crn/reaction1.crn --release-cutoff-1-1 4

class formal(s) = "d1a d1b d1c d2a d2b d2c d3a d3b d3c ? ? ?"
                | " .   .   .   .   .   .   .   .   .  ? ? ?"
    where {
        d1a = short() ; 
        d1b = short() ;
        d1c = short() ;
        d2a = short() ; 
        d2b = short() ;
        d2c = short() ;
        d3a = short() ; 
        d3b = short() ;
        d3c = short() };

module bifuel1(A,B,X,Y) = 
       "B2a B2b X2a X2b X2c X3a X3b X3c A2c* A2b* A2a* A1c* A1b* A1a* th* + th A1a A1b A1c A2a A2b Y2a Y2b Y2c Y3a Y3b + X3c* X3b* X3a* X2c* X2b* X2a* B2b* B2a* B1c*"
     | " (   (   (   (   (   (   (   (   .    (    (    (    (    (   (   + )   )   )   )   )   )   .   .   .   .   .  +  )    )    )    )    )    )    )    )    ."
    where {
        th = short();
    	A1a = A.d1a;
    	A1b = A.d1b;
    	A1c = A.d1c;
    	A2a = A.d2a;
    	A2b = A.d2b;
    	A2c = A.d2c;
	B1c = B.d1c;
	B2a = B.d2a;
	B2b = B.d2b;
	X2a = X.d2a;
	X2b = X.d2b;
	X2c = X.d2c;
	X3a = X.d3a;
	X3b = X.d3b;
	X3c = X.d3c;
	Y2a = Y.d2a;
	Y2b = Y.d2b;
	Y2c = Y.d2c;
	Y3a = Y.d3a;
	Y3b = Y.d3b };

module bifuel2(A,B,X,Y) = 
       "X1a X1b X1c X2a X2b X2c X3a X3b X3c A3a A3b A3c + A3c* A3b* A3a* X3c* X3b* X3a* X2c* X2b* X2a* B2b*" 
     | " .   .   .   (   (   (   (   (   (   (   (   (  +  )    )    )    )    )    )    )    )    )    ."
    where {
    	A3a = A.d3a;
    	A3b = A.d3b;
    	A3c = A.d3c;
	B2b = B.d2b;
	X1a = X.d1a;
	X1b = X.d1b;
	X1c = X.d1c;
	X2a = X.d2a;
	X2b = X.d2b;
	X2c = X.d2c;
	X3a = X.d3a;
	X3b = X.d3b;
	X3c = X.d3c };

module bimolecular(r) = infty(bifuel1(reactantA,reactantB,productX,productY)) 
       		      + infty(bifuel1(reactantB,reactantA,productY,productX)) 
		      + infty(bifuel2(reactantA,reactantB,productX,productY)) 
       		      + infty(bifuel2(reactantB,reactantA,productY,productX)) 
    where {
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
