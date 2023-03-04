=head1 NAME

dom_pretty_print.pm - functions for writing XML::DOM::Document object

=head1 SYNOPSIS

 usage: pretty_print($doc)
        $doc is XML::DOM::Document object
 or: traverse($doc, $filehandle)
     has to declare $docindex=-1; $lindex=0; first.

=head1 DESCRIPTION

  functions exist for the following
 
=head2 Methods

=over 12

=item C<pretty_print>

param
  $doc is XML::DOM::Document object

=item C<traverse>

param
  $doc is XML::DOM::Document object
  $filehandle is xml file
  calling program has to declare $docindex=-1; $lindex=0; first.

=back

=head1 AUTHOR

Dave Emmert - emmert@morgan.harvard.edu

=head1 SEE ALSO

WriteChado,  XML::DOM

=cut

#functions for writing XML::DOM::Document object to pretty XML 
#
# usage: pretty_print($doc)
#        $doc is XML::DOM::Document object
# or: traverse($doc, $filehandle)
#     has to declare $docindex=-1; $lindex=0; first.
#
%DecodeDefaultEntity =
(
 '"' => "&quot;",
 ">" => "&gt;",
 "<" => "&lt;",
 "'" => "&apos;",
 "&" => "&amp;"
);


sub pretty_print
{
  my ($node)=@_;
  $docindex=-1;
  $lindex=0;
  open(OUT, ">>test.xml");
  &traverse($node,\*OUT);
  close(OUT);
}
sub traverse {
  my($node,$filehandle)= @_;
  if ($node->getNodeType == ELEMENT_NODE) {
    $docindex++;
    if($lindex==0)
      { print $filehandle "\n";}
    my $attrs=$node->getAttributes();
    my @lists= $attrs->getValues;
    print  $filehandle ' 'x$docindex,"<",$node->getNodeName;
    foreach $attr(@lists)
      {
	print $filehandle ' ', $attr->getName,'=\'', $attr->getValue,'\'';
      }
    print $filehandle ">";

    foreach my $child ($node->getChildNodes()) {
      traverse($child,$filehandle);
    }
    if($lindex==0)
      {
	print $filehandle "\n",' 'x$docindex;
      }
    print $filehandle "</", $node->getNodeName, ">";
    $lindex=0;
    $docindex--;
  } elsif ($node->getNodeType() == TEXT_NODE) {
    print $filehandle &encodeText($node->getData, '<&>"');
    $lindex=1;
  }
  else {
    foreach my $child ($node->getChildNodes()){
     # print $node->getNodeType(),"\n";
      traverse($child,$filehandle);
  }
  }
}
sub encodeText
{
    my ($str, $default) = @_;
    return undef unless defined $str;

    if ($] >= 5.006) {
      $str =~ s/([$default])|(]]>)/
        defined ($1) ? $DecodeDefaultEntity{$1} : "]]&gt;" /egs;
    }
    else {
      $str =~ s/([\xC0-\xDF].|[\xE0-\xEF]..|[\xF0-\xFF]...)|([$default])|(]]>)/
        defined($1) ? XmlUtf8Decode ($1) :
        defined ($2) ? $DecodeDefaultEntity{$2} : "]]&gt;" /egs;
    }

#?? could there be references that should not be expanded?
# e.g. should not replace &#nn; &#xAF; and &abc;
#    $str =~ s/&(?!($ReName|#[0-9]+|#x[0-9a-fA-F]+);)/&amp;/go;

    $str;
}

sub XmlUtf8Decode
{
    my ($str, $hex) = @_;
    my $len = length ($str);
    my $n;

    if ($len == 2)
    {
        my @n = unpack "C2", $str;
        $n = (($n[0] & 0x3f) << 6) + ($n[1] & 0x3f);
    }
    elsif ($len == 3)
    {
        my @n = unpack "C3", $str;
        $n = (($n[0] & 0x1f) << 12) + (($n[1] & 0x3f) << 6) + 
                ($n[2] & 0x3f);
    }
    elsif ($len == 4)
    {
        my @n = unpack "C4", $str;
        $n = (($n[0] & 0x0f) << 18) + (($n[1] & 0x3f) << 12) + 
                (($n[2] & 0x3f) << 6) + ($n[3] & 0x3f);
    }
   elsif ($len == 1)   # just to be complete...
    {
        $n = ord ($str);
    }
    else
    {
        print "bad value [$str] for XmlUtf8Decode";
    }
    $hex ? sprintf ("&#x%x;", $n) : "&#$n;";
}

1;
