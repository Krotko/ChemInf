declare

  smi varchar2(512) := 'C1OC23COC45COC11COC67COC8(COC9(CO2)COC(CO1)(CO6)OCC(CO9)(OC4)OCC(CO5)(OC7)OC8)OC3';
  atoms Cheminf.atoms := Cheminf.atoms();
  bonds Cheminf.bonds := Cheminf.bonds();
  map matrix := matrix();
  ranks matrix := matrix();
  dist numarray := numarray();
  i binary_integer;
  j binary_integer;
  k binary_integer;
  part binary_integer;
  prevpart binary_integer;
  r binary_integer;
  tiedrank binary_integer;
  outp varchar2(1024);
  minsmi varchar2(512);

begin

  dbms_output.put_line (smi); dbms_output.put_line ('');

  smilin (smi, atoms, bonds);
  map := map_count (atoms, bonds);

  dbms_output.put_line ('Atoms:'); dbms_output.put_line ('');

  for i in 1..atoms.count loop

    outp := '';

    outp := outp || lpad (to_char (atoms(i).number), 3, ' ') || chr (9) || lpad (nvl (to_char (atoms(i).mass), ' '), 3, ' ') || chr (9) ||
    lpad (to_char (atoms(i).charge), 2, ' ');

    outp := outp || chr (9) || chr (9) || chr (9);

    outp := outp || to_char (atoms(i).degree) || chr (9) || to_char (atoms(i).hcount) || chr (9) || to_char (atoms(i).total_hcount) || chr (9)
    || to_char (atoms(i).connectivity) || chr (9) || to_char (atoms(i).valence);

    outp := outp || chr (9) || chr (9) || chr (9);

    outp := outp || to_char (atoms(i).ring_connectivity) || chr (9) || to_char (atoms(i).ring_membership) || chr (9) || lpad (nvl (to_char (atoms(i).ring_size), ' '), 2, ' ')
    || chr (9) || lpad (to_char (atoms(i).max_bonds_ringsize), 2, ' ') || chr (9) || to_char (atoms(i).aromaticity) || chr (9) || lpad (to_char (atoms(i).distance_inv), 6, ' ')
    || chr (9) || to_char (atoms(i).parity) || chr (9) || to_char (atoms(i).chirality) || chr (9) || lpad (to_char (atoms(i).symmetry_class), 2, ' ') || chr (9)
    || to_char (atoms(i).stereocenter) || chr (9) || lpad (to_char (atoms(i).rank), 2, ' ');

    dbms_output.put_line (outp);

  end loop;

  dbms_output.put_line (''); dbms_output.put_line ('Atoms count: '|| to_char (atoms.count) || '. Symmetry classes: ' || to_char (symmetry_classes_count (atoms))
  || '. Stereocenters: ' || to_char (stereocenters_count (atoms) || '. Chirality: ' || to_char (chiralatoms_count (atoms))));

  dbms_output.put_line (''); dbms_output.put_line ('Bonds:'); dbms_output.put_line ('');

  for i in 1..bonds.count loop

    outp := '';

    outp := outp || lpad (to_char (bonds(i).atom1), 2, ' ') || chr (9) || lpad (to_char (bonds(i).atom2), 2, ' ') || chr (9) || to_char (bonds(i).bondtype) ||
    chr (9) || chr (9) || chr (9) || to_char (bonds(i).ring_membership) || chr (9) || lpad (nvl (to_char (bonds(i).ring_size), ' '), 2, ' ') || chr (9) ||
    to_char (bonds(i).rotbond);

    dbms_output.put_line (outp);

  end loop;

  dbms_output.put_line (''); dbms_output.put_line ('Bonds count: '|| to_char (bonds.count) || '. Components: ' || to_char (components_count (atoms, bonds, map))
  || '. Rings: ' || to_char (rings_count (atoms, bonds)));

  dbms_output.put_line (''); dbms_output.put_line ('Canonical Codes:'); dbms_output.put_line ('');

  r := atoms(1).rank;

  for i in 2..atoms.count loop

    if atoms(i).rank > r then r := atoms(i).rank; end if;

  end loop;

  ranks.extend; ranks(1) := numarray();
  for i in 1..atoms.count loop ranks(1).extend; ranks(1)(i) := atoms(i).rank; end loop;
  dist.extend; dist(1) := r; part := 1; prevpart := 1;

  loop

    for i in prevpart..ranks.count loop

      if dist(i) < atoms.count then

         r := dist(i); tiedrank := 0;

         while tiedrank = 0 and r > 0 loop

           j := 1; k := 0;

           while k <= 1 and j <= atoms.count loop

             if ranks(i)(j) = r then k := k + 1; end if;

             if k = 2 and atoms(j).ring_connectivity > 0 then tiedrank := r; end if;

             j := j + 1;

           end loop;

           r := r - 1;

         end loop;

         for j in 1..atoms.count loop

           if ranks(i)(j) = tiedrank then

              ranks.extend; ranks(ranks.count) := numarray(); dist.extend;

              -- Breaking Ties
              for k in 1..atoms.count loop ranks(ranks.count).extend; ranks(ranks.count)(k) := 2 * ranks(i)(k); end loop;
              ranks(ranks.count)(j) := ranks(ranks.count)(j) - 1;

              for k in 1..atoms.count loop atoms(k).rank := ranks(ranks.count)(k); end loop;
              canon (atoms, bonds, map, dist(dist.count));
              for k in 1..atoms.count loop ranks(ranks.count)(k) := atoms(k).rank; end loop;

           end if;

         end loop;

      end if;

    end loop;

    if ranks.count = part then exit;
                          else prevpart := part + 1; part := ranks.count;
    end if;

  end loop;

  k := 0;

  minsmi := cangen (atoms, bonds);

  -- CCAP - Canonical Code Generation by Automorphism Permutation
  for i in prevpart..ranks.count loop

    for j in 1..atoms.count loop atoms(j).rank := ranks(i)(j); end loop;

    smi := cangen (atoms, bonds);

    dbms_output.put_line (smi);

    if smi < minsmi then minsmi := smi; end if;

  end loop;

  dbms_output.put_line (''); dbms_output.put_line ('Minimal Canonical Code: ' || minsmi || ' :' || to_char (ranks.count - prevpart + 1));

end;
