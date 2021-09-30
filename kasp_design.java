  private void PrimersDesign_KASP(int km, int i, int d, String CommnadsLines, int NpcrCombShow) {
        int ml = minprimerlength;   //  maximal length of additional 3'-ASP
        int xl = 0;                 //  calculated extended length
        int nv = 2;                 //  number SNP variants
        int l0 = S1seq[i].length(); //  seq length

        int z = nPRcombinations > 1 ? nPRcombinations : 1; // number sets to shown
        int mn = (km < ml) ? km : ml - 1;

        String p = "";
        String r = "";

        int x1 = S1oseq[i].indexOf('[');
        int x2 = S1oseq[i].indexOf(']', x1);

        if (x1 == -1 | x2 == -1) {
            return;
        }

        prminln = primerminlen;
        prmaxln = primermaxlen;
        prmintm = primermintm;
        prmaxtm = primermaxtm;
        pr5end = p5end;
        pr3end = p3end;
        dmr = 1;

        int CombShow = NpcrCombShow;
        if (CombShow < 1) {
            CombShow = 2;
        }

        ReadingControls(CommnadsLines, "p");
        String[] tg = new String[]{"", ""};    // splitted original sequence
        tg[0] = OnlyDNARead(S1oseq[i].substring(0, x1));
        tg[1] = OnlyDNARead(S1oseq[i].substring(x2));
        if (d == 2) {
            String s = tg[0];
            tg[0] = dna.ComplementDNA(tg[1]);
            tg[1] = dna.ComplementDNA(s);
        }
        int l1 = tg[0].length();
        int l2 = tg[1].length();
        int l3 = l0 - l1 - l2;

        if (l1 < 12 | l2 < 12) {
            return;
        }

        String[] q = pr5end.split("/");
        String[] s5 = new String[]{"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""};  // 5'tail     
        for (int g = 0; g < Math.min(q.length, s5.length); g++) {
            s5[g] = OnlyDNARead(q[g].toLowerCase());
        }

        String p51 = S1oseq[i].substring(x1 + 1, x2).toLowerCase();
        String[] st = new String[]{"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""}; // SNP   
        String[] t = p51.split("/");      // SNP 

        if (p51.contains("/")) {// [a/t]  [a/]
            nv = (int) p51.chars().filter(ch -> ch == '/').count() + 1;
        } else {                //  [s]
            int n = t[0].charAt(0);
            nv = t3[n][0];
            t = new String[nv];
            for (int g = 0; g < nv; g++) {
                t[g] = Character.toString((char) (t3[n][g + 1]));
            }
        }

        for (int h = 0; h < t.length; h++) {
            st[h] = OnlyDNARead(t[h]);
            if (d == 2) {
                st[h] = dna.ComplementDNA(t[h]);
            }
            if (st[h].length() < ml) {
                st[h] = t[h] + tg[1].substring(0, ml - t[h].length());
            }
        }
        for (int h = t.length; h < nv; h++) {
            if (st[h].length() < ml) {
                st[h] = tg[1].substring(0, ml);
            }
        }

        // SNP extending with target seq
        for (int x = 1; x < ml; x++) {
            int q0 = 0;
            for (int h = 0; h < nv - 1; h++) {
                p = st[h].substring(0, x);
                for (int n = h + 1; n < nv; n++) {
                    r = st[n].substring(0, x);
                    if (p == null ? r == null : p.equals(r)) {
                        h = nv;
                        q0 = 1;
                        break;
                    }
                }
            }
            if (q0 == 0) {
                xl = x + mn - 1;
                if (xl > ml) {
                    xl = ml;
                }
                break;
            }
        }

        if (xl > 0) {
            for (int n = 0; n < nv; n++) {
                st[n] = st[n].substring(0, xl);
            }
        }

        if (nv < 2 | xl == 0) {
            return;
        }

        //    SNP target Forward Primers
        PrimerData[] fpr = new PrimerData[nv];
        int fp = -1;
        int tm0 = (int) prmintm;
        for (int n = 0; n < nv; n++) {
            for (int j = 0; j < (1 + prmaxln - prminln); j++) {
                int x = l1 - prminln - j;
                p = tg[0].substring(x) + st[n];
                primer prm = new primer(p, 0);
                int tm = (int) prm.getTm();
                if (tm >= tm0) {
                    if (n == 0) {
                        tm0 = tm;
                    }
                    fpr[n] = new PrimerData();
                    fpr[n].x1t = x;
                    fpr[n].x2t = l1;
                    fpr[n].t = d;
                    fpr[n].tm = prm.getTm();
                    fpr[n].cg = prm.getGC();
                    fpr[n].dG = prm.getdG();
                    fpr[n].l = prm.getLC();
                    fpr[n].s = s5[n] + p;
                    fpr[n].ln = p.length();
                    fpr[n].q = Oligo.PrimerQuality(p);

                    if (d == 2) {
                        fpr[n].p = l2 + xl + fpr[n].ln;
                        fpr[n].n = "R_2_" + (fpr[n].p - fpr[n].ln + 1) + "-" + fpr[n].p;
                    } else {
                        fpr[n].p = 1 + l1 - fpr[n].ln;
                        fpr[n].n = "F_1_" + fpr[n].p + "-" + (fpr[n].p + fpr[n].ln - 1);
                    }
                    fp++;
                    if (dmr > 0) {
                        if (nPrimers > 0) {
                            if (Oligo.DimerLookList1(P1seq, fpr[n].s, min3dimer + dmr, -1) > 0) {
                                fp--;
                                break;
                            }
                        }
                        for (int h = 0; h <= n; h++) {
                            if (Oligo.DimerLook(fpr[n].s, fpr[h].s, min3dimer + dmr, -1) > 0) {
                                fp--;
                                break;
                            }
                        }
                    }
                    break;
                }
            }
            if (fp < n) {
                break;
            }
        }
        if (fp + 1 < nv) {
            return;
        }

//########################### SECOND REVERSE PRIMERS DESIGN  
        prminln = primerminlen;
        prmaxln = primermaxlen;
        prmintm = primermintm;
        prmaxtm = primermaxtm;
        provlp = OverlapPrimer;
        pr5end = p5end;
        pr3end = p3end;
        dmr = 1;
        ReadingControls(CommnadsLines, "r");
        PrimersList rplist = new PrimersList();
        rplist.setPrimerType(1);   //   1- Reverse
        RpT = new TaskCollector(); //   Reading Tasks

        if (d == 1) {
            RpT.insert(l0 - l2 + xl + l3, l0, 1, i, 0);
            RunPrimersDesign(false, rplist, RpT, i, "");
        }
        if (d == 2) {
            RpT.insert(0, l1 - xl - ml, 1, i, 0);
            RunPrimersDesign(true, rplist, RpT, i, "");
        }
        if (rplist.Amount() == 0) {
            return;
        }
//########################### END REVERSE PRIMERS DESIGN

// Combining primer to sets
        if (rplist.Amount() > 0) {
            PCRassayCollector pcrassay = new PCRassayCollector(null, null, i, "");
            pcrassay.PCRproductSize(CommnadsLines, l0);
            for (int j = 0; j < rplist.Amount(); j++) {
                if (CombShow < 1) {
                    break;
                }

                int x = 1;
                if (dmr > 0) {
                    for (int h = 0; h < nv; h++) {
                        if (Oligo.DimerLook(fpr[h].s, rplist.getPrimer(j), min3dimer + dmr, -1) > 0) {
                            x = 0;
                            break;
                        }
                    }
                }
                if (x == 1) {
                    int pcrl;
                    int ta;
                    if (d == 1) {
                        rplist.setLocation(j, 1 + l0 - rplist.getLocation(j));
                        rplist.setPrimerName(j, i, l0);
                        p = "R_1_" + (rplist.getLocation(j) - rplist.getLenght(j) + 1) + "-" + rplist.getLocation(j);     // p = rplist.getPrimerName(j);
                        pcrl = rplist.getLocation(j) + rplist.getLenght(j) - fpr[0].ln;
                        ta = (int) (Math.min(rplist.getTm(j), fpr[0].tm) + Math.log(pcrl));
                    } else {
                        rplist.setLocation(j, rplist.getLocation(j));
                        rplist.setPrimerName(j, i, l0);
                        p = "F_2_" + rplist.getLocation(j) + "-" + (rplist.getLocation(j) + rplist.getLenght(j) - 1);
                        pcrl = fpr[0].p - rplist.getLocation(j) + fpr[0].ln;
                        ta = (int) (Math.min(rplist.getTm(j), fpr[0].tm) + Math.log(pcrl));
                    }
                    if (pcrassay.getMaxPCRlen() >= pcrl) {
                        CombShow--;
                        z--;
                        for (int h = 0; h < nv; h++) {
                            sr.append(fpr[h].n).append("\t").append(fpr[h].s).append("\t").append(fpr[h].ln).append("\t").append(tools.NumberToSeq(fpr[h].tm, 1)).append("\t").append(tools.NumberToSeq(fpr[h].dG, 1)).append("\t").append(tools.NumberToSeq(fpr[h].cg, 1)).append("\t").append(fpr[h].l).append("\t").append(fpr[h].q).append("\n");
                        }
                        sr.append(p).append("\t").append(rplist.getPrimer(j)).append("\t").append(rplist.getLenght(j)).append("\t").append(tools.NumberToSeq(rplist.getTm(j), 1)).append("\t").append(tools.NumberToSeq(rplist.getdG(j), 1)).append("\t").append(tools.NumberToSeq(rplist.getCG(j), 1)).append("\t").append(rplist.getLinguisticComplexity(j)).append("\t").append(rplist.getPrimerQuality(j)).append("\t").append(pcrl).append("\t").append(ta).append("\n\n");
                        if (z < 0) {
                            break;
                        }
                    }
                }
            }
        }

    }
