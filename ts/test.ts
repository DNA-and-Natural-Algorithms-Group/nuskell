function formal(s) =
    "? a b c" | "? . . ."
    where {
        a = short() ;
        b = long() ;
        c = short()
    };
function Bi(x1, x2) =
    "d2 d3 d4" |
    " . . ."
    where {
        d2 = x1.b ;
        d3 = x1.c ;
        d4 = x2.a
    };
function main(crn) = infty(Bi(crn[0].reactants[0], crn[0].reactants[1]))
