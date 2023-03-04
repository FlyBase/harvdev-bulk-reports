sub conversions {

sub convers {
# Does necessary conversions on all input strings

local($string) = shift;


	chop($string);

	$string =~ s/\'/\'\'/g;
	$string =~ s/^\s+|^\n|\n$|\n!$|\n\s$//g;

    $string =~ s/\$a/&agr\;/g;
    $string =~ s/\$A/&Agr\;/g;
    $string =~ s/\$b/&bgr\;/g;
    $string =~ s/\$B/&Bgr\;/g;
    $string =~ s/\$g/&ggr\;/g;
    $string =~ s/\$G/&Ggr\;/g;
    $string =~ s/\$d/&dgr\;/g;
    $string =~ s/\$D/&Dgr\;/g;
    $string =~ s/\$ee/&eegr\;/g;
    $string =~ s/\$EE/&EEgr\;/g;
    $string =~ s/\$e/&egr\;/g;
    $string =~ s/\$E/&Egr\;/g;
    $string =~ s/\$i/&igr\;/g;
    $string =~ s/\$I/&Igr\;/g;
    $string =~ s/\$kh/&khgr\;/g;
    $string =~ s/\$KH/&KHgr\;/g;
    $string =~ s/\$k/&kgr\;/g;
    $string =~ s/\$K/&Kgr\;/g;
    $string =~ s/\$l/&lgr\;/g;
    $string =~ s/\$L/&Lgr\;/g;
    $string =~ s/\$m/&mgr\;/g;
    $string =~ s/\$M/&Mgr\;/g;
    $string =~ s/\$n/&ngr\;/g;
    $string =~ s/\$N/&Ngr\;/g;
    $string =~ s/\$oh/&ohgr\;/g;
    $string =~ s/\$OH/&OHgr\;/g;
    $string =~ s/\$o/&ogr\;/g;
    $string =~ s/\$O/&Ogr\;/g;
    $string =~ s/\$ph/&phgr\;/g;
    $string =~ s/\$PH/&PHgr\;/g;
    $string =~ s/\$ps/&psgr\;/g;
    $string =~ s/\$PS/&PSgr\;/g;
    $string =~ s/\$p/&pgr\;/g;
    $string =~ s/\$P/&Pgr\;/g;
    $string =~ s/\$r/&rgr\;/g;
    $string =~ s/\$R/&Rgr\;/g;
    $string =~ s/\$s/&sgr\;/g;
    $string =~ s/\$S/&Sgr\;/g;
    $string =~ s/\$th/&thgr\;/g;
    $string =~ s/\$TH/&THgr\;/g;
    $string =~ s/\$t/&tgr\;/g;
    $string =~ s/\$T/&Tgr\;/g;
    $string =~ s/\$u/&ugr\;/g;
    $string =~ s/\$U/&Ugr\;/g;
    $string =~ s/\$z/&zgr\;/g;
    $string =~ s/\$Z/&Zgr\;/g;
    $string =~ s/\$x/&xgr\;/g;
    $string =~ s/\$X/&Xgr\;/g;
    $string =~ s/\]{2,2}/\<\/down\>/g;
    $string =~ s/\[{2,2}/\<down\>/g;
    $string =~ s/\[/\<up\>/g;
    $string =~ s/\]/\<\/up\>/g;

    return($string);
}

sub con_greeks {
# Converts greek characters in symbols

    local($string) = shift;


    $string =~ s/\$a/&agr\;/g;
    $string =~ s/\$A/&Agr\;/g;
    $string =~ s/\$b/&bgr\;/g;
    $string =~ s/\$B/&Bgr\;/g;
    $string =~ s/\$g/&ggr\;/g;
    $string =~ s/\$G/&Ggr\;/g;
    $string =~ s/\$d/&dgr\;/g;
    $string =~ s/\$D/&Dgr\;/g;
    $string =~ s/\$ee/&eegr\;/g;
    $string =~ s/\$EE/&EEgr\;/g;
    $string =~ s/\$e/&egr\;/g;
    $string =~ s/\$E/&Egr\;/g;
    $string =~ s/\$i/&igr\;/g;
    $string =~ s/\$I/&Igr\;/g;
    $string =~ s/\$th/&thgr\;/g;
    $string =~ s/\$TH/&THgr\;/g;
    $string =~ s/\$k/&kgr\;/g;
    $string =~ s/\$K/&Kgr\;/g;
    $string =~ s/\$l/&lgr\;/g;
    $string =~ s/\$L/&Lgr\;/g;
    $string =~ s/\$m/&mgr\;/g;
    $string =~ s/\$M/&Mgr\;/g;
    $string =~ s/\$n/&ngr\;/g;
    $string =~ s/\$N/&Ngr\;/g;
    $string =~ s/\$oh/&ohgr\;/g;
    $string =~ s/\$OH/&OHgr\;/g;
    $string =~ s/\$o/&ogr\;/g;
    $string =~ s/\$O/&Ogr\;/g;
    $string =~ s/\$ph/&phgr\;/g;
    $string =~ s/\$PH/&PHgr\;/g;
    $string =~ s/\$ps/&psgr\;/g;
    $string =~ s/\$PS/&PSgr\;/g;
    $string =~ s/\$p/&pgr\;/g;
    $string =~ s/\$P/&Pgr\;/g;
    $string =~ s/\$r/&rgr\;/g;
    $string =~ s/\$R/&Rgr\;/g;
    $string =~ s/\$s/&sgr\;/g;
    $string =~ s/\$S/&Sgr\;/g;
    $string =~ s/\$th/&thgr\;/g;
    $string =~ s/\$TH/&THgr\;/g;
    $string =~ s/\$t/&tgr\;/g;
    $string =~ s/\$T/&Tgr\;/g;
    $string =~ s/\$u/&ugr\;/g;
    $string =~ s/\$U/&Ugr\;/g;
    $string =~ s/\$z/&zgr\;/g;
    $string =~ s/\$Z/&Zgr\;/g;
    $string =~ s/\$x/&xgr\;/g;
    $string =~ s/\$X/&Xgr\;/g;
    $string =~ s/\]{2,2}/\<\/down\>/g;
    $string =~ s/\[{2,2}/\<down\>/g;
    $string =~ s/\[/\<up\>/g;
    $string =~ s/\]/\<\/up\>/g;
    
    return($string);
}

sub decon {
# Converts SGML-formatted symbols to 'symbol_plain' format (modified from conv_greeks)

    local($string) = shift;

    $string =~ s/&agr\;/alpha/g;
    $string =~ s/&Agr\;/Alpha/g;
    $string =~ s/&bgr\;/beta/g;
    $string =~ s/&Bgr\;/Beta/g;
    $string =~ s/&ggr\;/gamma/g;
    $string =~ s/&Ggr\;/Gamma/g;
    $string =~ s/&dgr\;/delta/g;
    $string =~ s/&Dgr\;/Delta/g;
    $string =~ s/&egr\;/epsilon/g;
    $string =~ s/&Egr\;/Epsilon/g;
    $string =~ s/&zgr\;/zeta/g;
    $string =~ s/&Zgr\;/Zeta/g;
    $string =~ s/&eegr\;/eta/g;
    $string =~ s/&EEgr\;/Eta/g;
    $string =~ s/&thgr\;/theta/g;
    $string =~ s/&THgr\;/Theta/g;
    $string =~ s/&igr\;/iota/g;
    $string =~ s/&Igr\;/Iota/g;
    $string =~ s/&kgr\;/kappa/g;
    $string =~ s/&Kgr\;/Kappa/g;
    $string =~ s/&lgr\;/lambda/g;
    $string =~ s/&Lgr\;/Lambda/g;
    $string =~ s/&mgr\;/mu/g;
    $string =~ s/&Mgr\;/Mu/g;
    $string =~ s/&ngr\;/nu/g;
    $string =~ s/&Ngr\;/Nu/g;
    $string =~ s/&xgr\;/xi/g;
    $string =~ s/&Xgr\;/Xi/g;
    $string =~ s/&ogr\;/omicron/g;
    $string =~ s/&Ogr\;/Omicron/g;
    $string =~ s/&pgr\;/pi/g;
    $string =~ s/&Pgr\;/Pi/g;
    $string =~ s/&rgr\;/rho/g;
    $string =~ s/&Rgr\;/Rho/g;
    $string =~ s/&sgr\;/sigma/g;
    $string =~ s/&Sgr\;/Sigma/g;
    $string =~ s/&tgr\;/tau/g;
    $string =~ s/&Tgr\;/Tau/g;
    $string =~ s/&ugr\;/upsilon/g;
    $string =~ s/&Ugr\;/Upsilon/g;
    $string =~ s/&phgr\;/phi/g;
    $string =~ s/&PHgr\;/Phi/g;
    $string =~ s/&khgr\;/chi/g;
    $string =~ s/&KHgr\;/Chi/g;
    $string =~ s/&psgr\;/psi/g;
    $string =~ s/&PSgr\;/Psi/g;
    $string =~ s/&ohgr\;/omega/g;
    $string =~ s/&OHgr\;/Omega/g;
    $string =~ s/\<\/down\>/\]\]/g;
    $string =~ s/\<down\>/\[\[/g;
    $string =~ s/\<up\>/\[/g;
    $string =~ s/\<\/up\>/\]/g;

    return($string);
}

	
sub recon {
# Converts symbol_plain symbols to SGML format (modified from conv_greeks)

    local($string) = shift;

    $string =~ s/alpha/&agr\;/g;
    $string =~ s/Alpha/&Agr\;/g;
    $string =~ s/beta/&bgr\;/g;
    $string =~ s/Beta/&Bgr\;/g;
    $string =~ s/gamma/&ggr\;/g;
    $string =~ s/Gamma/&Ggr\;/g;
    $string =~ s/delta/&dgr\;/g;
    $string =~ s/Delta/&Dgr\;/g;
    $string =~ s/epsilon/&egr\;/g;
    $string =~ s/Epsilon/&Egr\;/g;
    $string =~ s/zeta/&zgr\;/g;
    $string =~ s/Zeta/&Zgr\;/g;
    $string =~ s/eta/&eegr\;/g;
    $string =~ s/Eta/&EEgr\;/g;
    $string =~ s/theta/&thgr\;/g;
    $string =~ s/Theta/&THgr\;/g;
    $string =~ s/iota/&igr\;/g;
    $string =~ s/Iota/&Igr\;/g;
    $string =~ s/kappa/&kgr\;/g;
    $string =~ s/Kappa/&Kgr\;/g;
    $string =~ s/lambda/&lgr\;/g;
    $string =~ s/Lambda/&Lgr\;/g;
    $string =~ s/mu/&mgr\;/g;
    $string =~ s/Mu/&Mgr\;/g;
    $string =~ s/nu/&ngr\;/g;
    $string =~ s/Nu/&Ngr\;/g;
    $string =~ s/xi/&xgr\;/g;
    $string =~ s/Xi/&Xgr\;/g;
    $string =~ s/omicron/&ogr\;/g;
    $string =~ s/Omicron/&Ogr\;/g;
    $string =~ s/pi/&pgr\;/g;
    $string =~ s/Pi/&Pgr\;/g;
    $string =~ s/rho/&rgr\;/g;
    $string =~ s/Rho/&Rgr\;/g;
    $string =~ s/sigma/&sgr\;/g;
    $string =~ s/Sigma/&Sgr\;/g;
    $string =~ s/tau/&tgr\;/g;
    $string =~ s/Tau/&Tgr\;/g;
    $string =~ s/upsilon/&ugr\;/g;
    $string =~ s/Upsilon/&Ugr\;/g;
    $string =~ s/phi/&phgr\;/g;
    $string =~ s/Phi/&PHgr\;/g;
    $string =~ s/chi/&khgr\;/g;
    $string =~ s/Chi/&KHgr\;/g;
    $string =~ s/psi/&psgr\;/g;
    $string =~ s/Psi/&PSgr\;/g;
    $string =~ s/omega/&ohgr\;/g;
    $string =~ s/Omega/&OHgr\;/g;
    $string =~ s/\]\]/\<\/down\>/g;
    $string =~ s/\[\[/\<down\>/g;
    $string =~ s/\[/\<up\>/g;
    $string =~ s/\]/\<\/up\>/g;

    return($string);
}

sub htmlize {
## convert controlled characters to html-legal

    local($string) = shift;

    $string =~ s/\&/\&amp\;/g;

    $string =~ s/\</\&lt\;/g;
    $string =~ s/\>/\&gt\;/g;

    return($string);

}

sub utf2sgml {
    my ($string) = shift;

    $string =~ s/\x{03B1}/&agr\;/g;
    $string =~ s/\x{0391}/&Agr\;/g;
    $string =~ s/\x{03B2}/&bgr\;/g;
    $string =~ s/\x{0392}/&Bgr\;/g;
    $string =~ s/\x{03B3}/&ggr\;/g;
    $string =~ s/\x{0393}/&Ggr\;/g;
    $string =~ s/\x{03B4}/&dgr\;/g;
    $string =~ s/\x{0394}/&Dgr\;/g;
    $string =~ s/\x{03B5}/&egr\;/g;
    $string =~ s/\x{0395}/&Egr\;/g;
    $string =~ s/\x{03B6}/&zgr\;/g;
    $string =~ s/\x{0396}/&Zgr\;/g;
    $string =~ s/\x{03B7}/&eegr\;/g;
    $string =~ s/\x{0397}/&EEgr\;/g;
    $string =~ s/\x{03B8}/&thgr\;/g;
    $string =~ s/\x{0398}/&THgr\;/g;
    $string =~ s/\x{03B9}/&igr\;/g;
    $string =~ s/\x{0399}/&Igr\;/g;
    $string =~ s/\x{03BA}/&kgr\;/g;
    $string =~ s/\x{039A}/&Kgr\;/g;
    $string =~ s/\x{03BB}/&lgr\;/g;
    $string =~ s/\x{039B}/&Lgr\;/g;
    $string =~ s/\x{03BC}/&mgr\;/g;
    $string =~ s/\x{039C}/&Mgr\;/g;
    $string =~ s/\x{03BD}/&ngr\;/g;
    $string =~ s/\x{039D}/&Ngr\;/g;
    $string =~ s/\x{03BE}/&xgr\;/g;
    $string =~ s/\x{039E}/&Xgr\;/g;
    $string =~ s/\x{03BF}/&ogr\;/g;
    $string =~ s/\x{039F}/&Ogr\;/g;
    $string =~ s/\x{03C0}/&pgr\;/g;
    $string =~ s/\x{03A0}/&Pgr\;/g;
    $string =~ s/\x{03C1}/&rgr\;/g;
    $string =~ s/\x{03A1}/&Rgr\;/g;
    $string =~ s/\x{03C3}/&sgr\;/g;
    $string =~ s/\x{03A3}/&Sgr\;/g;
    $string =~ s/\x{03C4}/&tgr\;/g;
    $string =~ s/\x{03A4}/&Tgr\;/g;
    $string =~ s/\x{03C5}/&ugr\;/g;
    $string =~ s/\x{03A5}/&Ugr\;/g;
    $string =~ s/\x{03C6}/&phgr\;/g;
    $string =~ s/\x{03A6}/&PHgr\;/g;
    $string =~ s/\x{03C7}/&khgr\;/g;
    $string =~ s/\x{03A7}/&KHgr\;/g;
    $string =~ s/\x{03C8}/&psgr\;/g;
    $string =~ s/\x{03A8}/&PSgr\;/g;
    $string =~ s/\x{03C9}/&ohgr\;/g;
    $string =~ s/\x{03A9}/&OHgr\;/g;

    $string =~ s/\<\/down\>/\]\]/g;
    $string =~ s/\<down\>/\[\[/g;
    $string =~ s/\<up\>/\[/g;
    $string =~ s/\<\/up\>/\]/g;

    return ($string);

}

sub atp_utf2sgml {
    my ($string) = shift;

    $string =~ s/\x{03B1}/&agr\;/g;
    $string =~ s/\x{0391}/&Agr\;/g;
    $string =~ s/\x{03B2}/&bgr\;/g;
    $string =~ s/\x{0392}/&Bgr\;/g;
    $string =~ s/\x{03B3}/&ggr\;/g;
    $string =~ s/\x{0393}/&Ggr\;/g;
    $string =~ s/\x{03B4}/&dgr\;/g;
    $string =~ s/\x{0394}/&Dgr\;/g;
    $string =~ s/\x{03B5}/&egr\;/g;
    $string =~ s/\x{0395}/&Egr\;/g;
    $string =~ s/\x{03B6}/&zgr\;/g;
    $string =~ s/\x{0396}/&Zgr\;/g;
    $string =~ s/\x{03B7}/&eegr\;/g;
    $string =~ s/\x{0397}/&EEgr\;/g;
    $string =~ s/\x{03B8}/&thgr\;/g;
    $string =~ s/\x{0398}/&THgr\;/g;
    $string =~ s/\x{03B9}/&igr\;/g;
    $string =~ s/\x{0399}/&Igr\;/g;
    $string =~ s/\x{03BA}/&kgr\;/g;
    $string =~ s/\x{039A}/&Kgr\;/g;
    $string =~ s/\x{03BB}/&lgr\;/g;
    $string =~ s/\x{039B}/&Lgr\;/g;
    $string =~ s/\x{03BC}/&mgr\;/g;
    $string =~ s/\x{039C}/&Mgr\;/g;
    $string =~ s/\x{03BD}/&ngr\;/g;
    $string =~ s/\x{039D}/&Ngr\;/g;
    $string =~ s/\x{03BE}/&xgr\;/g;
    $string =~ s/\x{039E}/&Xgr\;/g;
    $string =~ s/\x{03BF}/&ogr\;/g;
    $string =~ s/\x{039F}/&Ogr\;/g;
    $string =~ s/\x{03C0}/&pgr\;/g;
    $string =~ s/\x{03A0}/&Pgr\;/g;
    $string =~ s/\x{03C1}/&rgr\;/g;
    $string =~ s/\x{03A1}/&Rgr\;/g;
    $string =~ s/\x{03C3}/&sgr\;/g;
    $string =~ s/\x{03A3}/&Sgr\;/g;
    $string =~ s/\x{03C4}/&tgr\;/g;
    $string =~ s/\x{03A4}/&Tgr\;/g;
    $string =~ s/\x{03C5}/&ugr\;/g;
    $string =~ s/\x{03A5}/&Ugr\;/g;
    $string =~ s/\x{03C6}/&phgr\;/g;
    $string =~ s/\x{03A6}/&PHgr\;/g;
    $string =~ s/\x{03C7}/&khgr\;/g;
    $string =~ s/\x{03A7}/&KHgr\;/g;
    $string =~ s/\x{03C8}/&psgr\;/g;
    $string =~ s/\x{03A8}/&PSgr\;/g;
    $string =~ s/\x{03C9}/&ohgr\;/g;
    $string =~ s/\x{03A9}/&OHgr\;/g;

    return ($string);

}

sub trim {
  my @out = @_;
  for (@out) {
    s/^\s+//;
    s/\s+$//;
  }
  return wantarray ? @out : $out[0];
}

  
=head2 spell_greek

takes a word as a parameter and spells out any greek symbols encoded
within (eg s/&agr;/alpha/g)

=cut

sub spell_greek
{
    my $name = shift;

    $name =~ s/&agr\;/alpha/g;
    $name =~ s/&Agr\;/Alpha/g;
    $name =~ s/&bgr\;/beta/g;
    $name =~ s/&Bgr\;/Beta/g;
    $name =~ s/&ggr\;/gamma/g;
    $name =~ s/&Ggr\;/Gamma/g;
    $name =~ s/&dgr\;/delta/g;
    $name =~ s/&Dgr\;/Delta/g;
    $name =~ s/&egr\;/epsilon/g;
    $name =~ s/&Egr\;/Epsilon/g;
    $name =~ s/&zgr\;/zeta/g;
    $name =~ s/&Zgr\;/Zeta/g;
    $name =~ s/&eegr\;/eta/g;
    $name =~ s/&EEgr\;/Eta/g;
    $name =~ s/&thgr\;/theta/g;
    $name =~ s/&THgr\;/Theta/g;
    $name =~ s/&igr\;/iota/g;
    $name =~ s/&Igr\;/Iota/g;
    $name =~ s/&kgr\;/kappa/g;
    $name =~ s/&Kgr\;/Kappa/g;
    $name =~ s/&lgr\;/lambda/g;
    $name =~ s/&Lgr\;/Lambda/g;
    $name =~ s/&mgr\;/mu/g;
    $name =~ s/&Mgr\;/Mu/g;
    $name =~ s/&ngr\;/nu/g;
    $name =~ s/&Ngr\;/Nu/g;
    $name =~ s/&xgr\;/xi/g;
    $name =~ s/&Xgr\;/Xi/g;
    $name =~ s/&ogr\;/omicron/g;
    $name =~ s/&Ogr\;/Omicron/g;
    $name =~ s/&pgr\;/pi/g;
    $name =~ s/&Pgr\;/Pi/g;
    $name =~ s/&rgr\;/rho/g;
    $name =~ s/&Rgr\;/Rho/g;
    $name =~ s/&sgr\;/sigma/g;
    $name =~ s/&Sgr\;/Sigma/g;
    $name =~ s/&tgr\;/tau/g;
    $name =~ s/&Tgr\;/Tau/g;
    $name =~ s/&ugr\;/upsilon/g;
    $name =~ s/&Ugr\;/Upsilon/g;
    $name =~ s/&phgr\;/phi/g;
    $name =~ s/&PHgr\;/Phi/g;
    $name =~ s/&khgr\;/chi/g;
    $name =~ s/&KHgr\;/Chi/g;
    $name =~ s/&psgr\;/psi/g;
    $name =~ s/&PSgr\;/Psi/g;
    $name =~ s/&ohgr\;/omega/g;
    $name =~ s/&OHgr\;/Omega/g;
    $name =~ s/<up>/\[/g;
    $name =~ s/<\/up>/\]/g;
    $name =~ s/<down>/\[\[/g;
    $name =~ s/<\/down>/\]\]/g;

    return $name;
}


}

1;
