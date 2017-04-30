#
# The section 6, figure 9 "xchange gates" from Thachuk, Winfree, Soloveichik, LNCS xxx, 2015
# 
# Coded by Erik Winfree (winfree@caltech.edu)
#
# _WF stands for "waste fuel".  There are two "wastes" produced by the standard reaction pathway that aren't formally waste, because they are not inert.
# So make them fuels, formally.  They are build by waste1().
# Unfortunately, although it appears "correct" in some sense, this scheme cannot be proved correct by bisimulation or pathway decomposition, even on just A+B->X+Y.
#
# To avoid combinatorial explosion, this must be invoked along the lines of:
# ./verify ts/thachuk_xchange_WF.ts crn/reaction1.crn --release-cutoff-1-1 4

class formal(s) = "d1a d1b d1c d2a d2b d2c d3a d3b d3c ? ? ?"
                | " .   .   .   .   .   .   .   .   .  . . ."
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

module bifuel1(A,B,X,Y,th) = 
       "B2a B2b X2a X2b X2c X3a X3b X3c A2c* A2b* A2a* A1c* A1b* A1a* th* + th A1a A1b A1c A2a A2b Y2a Y2b Y2c Y3a Y3b + X3c* X3b* X3a* X2c* X2b* X2a* B2b* B2a* B1c*"
     | " (   (   (   (   (   (   (   (   .    (    (    (    (    (   (   + )   )   )   )   )   )   .   .   .   .   .  +  )    )    )    )    )    )    )    )    ."
    where {
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

module bifuel2(A,B,X,Y,th) = 
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

module waste1(A,B,X,Y,th) = 
       "th B1a B1b B1c B2a B2b X2a X2b X2c X3a X3b + X3c* X3b* X3a* X2c* X2b* X2a* B2b* B2a* B1c*"
     | ".   .   .   (   (   (   (   (   (   (   (  +  .    )    )    )    )    )    )    )    )"
    where {
    	B1a = B.d1a;
    	B1b = B.d1b;
    	B1c = B.d1c;
	B2a = B.d2a;
	B2b = B.d2b;
	X2a = X.d2a;
	X2b = X.d2b;
	X2c = X.d2c;
	X3a = X.d3a;
	X3b = X.d3b;
	X3c = X.d3c };
        

module bimolecular(r) = infty(bifuel1(reactantA,reactantB,productX,productY,th)) 
       		      + infty(bifuel1(reactantB,reactantA,productY,productX,th)) 
		      + infty(bifuel2(reactantA,reactantB,productX,productY,th)) 
       		      + infty(bifuel2(reactantB,reactantA,productY,productX,th)) 
		      + infty(waste1(reactantA,reactantB,productX,productY,th)) 
       		      + infty(waste1(reactantB,reactantA,productY,productX,th)) 
    where {
        reactantA = r.reactants[0];
	reactantB = r.reactants[1];
	productX = r.products[0];
	productY = r.products[1];
	th = short() };

module rxn(r) =
    if len(r.reactants) == 2 and len(r.products) == 2 then
        bimolecular(r)
    else r[0];
# this last line is a hack to cause a crash when the target CRN is invalid for this scheme

module main(crn) = sum(map(rxn, crn))
    where
        crn = irrev_reactions(crn)
