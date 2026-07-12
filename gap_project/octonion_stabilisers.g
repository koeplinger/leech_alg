#############################################################################
##
##  octonion_stabilisers.g -- identify the three octonion-automorphism
##  lattice stabilisers used by the Aut(Lambda,+,star) computation.
##
##  Reads 8x8 rational generators produced by
##      python_project/src/verify_aut_octonion_crosscheck.py
##  and identifies the three finite subgroups of the compact G_2 = Aut(O, .s):
##
##      StabL      = { A : A(L)     = L     }   expected order  1344
##      StabLs     = { A : A(Ls)    = Ls    }   expected order   168
##      StabLsbar  = { A : A(Lsbar) = Lsbar }   expected order 12096
##
##  The last one is |G_2(2)|; the paper's Ls_bar, not L, is the lattice whose
##  octonion-automorphism stabiliser is the full G_2(2).
##
#############################################################################

Read("octonion_stabilisers_gens.g");

Ident := function(name, gens)
  local G;
  G := Group(gens);
  Print(name, ": order ", Size(G),
        ",  structure ", StructureDescription(G), "\n");
  if Size(G) = 12096 then
    Print("   isomorphic to G_2(2) = Aut(G_2(2)') : ",
          IsomorphismGroups(G, AtlasGroup("U3(3).2")) <> fail, "\n");
  fi;
  return G;
end;

G1 := Ident("StabL     ", StabL);
G2 := Ident("StabLs    ", StabLs);
G3 := Ident("StabLsbar ", StabLsbar);

Print("StabLs  <= StabLsbar : ", IsSubset(G3, G2), "\n");
Print("|StabL ^ StabLs ^ StabLsbar| = ",
      Size(Intersection(Intersection(G1,G2),G3)),
      "   (the blockwise automorphisms of Lambda)\n");
Print("structure of that intersection: ",
      StructureDescription(Intersection(Intersection(G1,G2),G3)), "\n");
