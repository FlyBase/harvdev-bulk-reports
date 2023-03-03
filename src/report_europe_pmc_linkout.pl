use strict;
use XML::DOM;
use DBI;
use Digest::MD5  qw(md5 md5_hex md5_base64);

=head
here for old ticket WEB_455, to generate linkout file for PMC so they can point to FlyBase for all FlyBase's pubmed record. query to retrieve all no obsoleted pubmed:
select d.accession from dbxref d, pub p, pub_dbxref pd, db d1 where p.pub_id=pd.pub_id and pd.dbxref_id=d.dbxref_id and d.db_id=d1.db_id and d1.name='pubmed' and p.is_obsolete='false';
NOTE - this script has been adapted to run in docker by GoCD.

generate two files to PMC:
1. profile.xml
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<providers>
  <provider>
    <id>1452</id>
    <resourceName>FlyBase </resourceName>
    <description>A Database of Drosophila Genes &amp; Genomes</description>
    <email>helpfb@morgan.harvard.edu</email>
  </provider>
</providers>

2. epmc-flybase.xml
<?xml version="1.0" encoding="UTF-8"?>
<links>
  <link providerId="1452">
    <resource>
      <url>http://flybase.org/cgi-bin/uniq.html?db=fbrf&amp;field=pubmed_id&amp;context=15926476</url>
    </resource>
    <record>
      <source>MED</source>
      <id>15926476</id>
    </record>
  </link>
  <link providerId="1452">
    <resource>
      <url>http://flybase.org/cgi-bin/uniq.html?db=fbrf&amp;field=pubmed_id&amp;context=15926510</url>
    </resource>
    <record>
      <source>MED</source>
      <id>15926510</id>
    </record>
  </link>
 ...
</links>

=cut

if (@ARGV != 4) {
    print "\n USAGE: $0 pg_server db_name pg_username pg_password\n\n";
    print "\teg: $0 flysql19 production_chado zhou pwd \n\n";
    exit;
}

my $server = shift(@ARGV);
my $db = shift(@ARGV);
my $user = shift(@ARGV);
my $pwd = shift(@ARGV);


#my $file_chadoxml=shift(@ARGV);
#default file name is epmc-flybase.xml
my $file_chadoxml='/src/output/epmc-flybase.xml';


if (-e $file_chadoxml) {
  unlink($file_chadoxml) or die "$file_chadoxml: $!"
}

my $file_chadoxml_gz= $file_chadoxml.".gz";
if (-e $file_chadoxml_gz) {
  unlink($file_chadoxml_gz) or die "$file_chadoxml_gz: $!"
}

open (OUT, ">$file_chadoxml") or die "unable to write to $file_chadoxml";
   print  OUT '<?xml version="1.0" encoding="UTF-8"?>'."\n<links>";

my $doc=new XML::DOM::Document();
my $header='<?xml version="1.0" encoding="UTF-8"?>';

my $dsource = sprintf("dbi:Pg:dbname=%s;host=%s;port=5432",$db,$server);
my $dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";


my $sql_FBrf1=sprintf("select d.accession from dbxref d, pub p, pub_dbxref pd, db d1 where p.pub_id=pd.pub_id and pd.dbxref_id=d.dbxref_id and d.db_id=d1.db_id and d1.name='pubmed' and p.is_obsolete='false' ");
#print "\n$sql_FBrf1";exit; 

my $q_FBrf1 = $dbh->prepare  ($sql_FBrf1);
$q_FBrf1->execute or die" CAN'T GET UNIPROT ACCESSIONS FROM CHADO:\n$sql_FBrf1\n";
 my ($accession );
my $providerId=1452;
my $source='MED';
my $url_prefix='http://flybase.org/cgi-bin/uniq.html?db=fbrf&amp;field=pubmed_id&amp;context=';

 while (($accession) = $q_FBrf1->fetchrow_array()) {    
     my $url=$url_prefix.$accession;
     my $node_link=&_get_link_node($accession, $source, $url, $providerId);
     if (defined $node_link && $node_link ne ""){
	 &_traverse($node_link);
     }
}

print OUT "\n</links>";
close (OUT);
system ("/bin/gzip -f $file_chadoxml");


# This perl ftp upload works on flysql servers.
# However, when used in GoCD, something gets corrupted.
# i.e., I get "Invalid compressed data--format violated" error when I try to 
#     download/gunzip the just-uploaded file from the EBI site.
# So, letting go-agent handle ftp transfer (with .netrc file).
################################################################################
# #her to upload into PMC server
# use warnings;
# use Net::FTP;

# my ($ftp, $host, $user, $pass, $dir, $fpath);

# $host = "labslink.ebi.ac.uk";
# $user = "elinks";
# $pass = "8VhrURVH";
# $dir = "c5odmhs7";

# $fpath = $file_chadoxml_gz;

# $ftp = Net::FTP->new($host, Debug => 1, Passive => 1) or die "Couldn't connect: $@\n";
# $ftp->login($user, $pass) || die  "problem1 here:". $ftp->message;
# $ftp->cwd($dir);
# if ( $ftp->delete ( $fpath ) ) {
#    print "Old File <$fpath> was deleted\n";
#    } else {
#    print "Old File <$fpath> failed to delete\n";
# }
# $ftp->put($fpath) || die "problem2 here:$fpath ".$ftp->message;
# $ftp->quit;

# print $ftp->message;
################################################################################


=header
create the link node for PMC file
=cut
sub _get_link_node(){
   my ($id, $source, $url, $providerId)=@_;
   if ($id eq "" || $source eq "" || $url eq "" || $providerId eq ""){
       return;
   }

   my $node_1=$doc->createTextNode($url);
   my $node_2=$doc->createTextNode($source);
   my $node_3=$doc->createTextNode($id);

   my $node_url=$doc->createElement('url');
   my $node_link=$doc->createElement('link');
   my $node_resource=$doc->createElement('resource');
   my $node_record=$doc->createElement('record');
   my $node_source=$doc->createElement('source');
   my $node_id=$doc->createElement('id');
   $node_url->appendChild($node_1);
   $node_source->appendChild($node_2);
   $node_id->appendChild($node_3);
   $node_resource->appendChild($node_url);
   $node_record->appendChild($node_source);
   $node_record->appendChild($node_id);
   $node_link->appendChild($node_resource);
   $node_link->appendChild($node_record);
   $node_link->setAttribute('providerId', $providerId);
   return $node_link;
};


=header 
   private method to print out the generated chadoxml
=cut
sub _traverse {
    my($node, $indent)= @_;
    if (!(defined $indent)){
	$indent=1;
    }
    if ($node->getNodeType == ELEMENT_NODE) {

      my $att_id=$node->getAttribute ('id');
      #$att_id=&_decon_xml($att_id);
      my $att_op=$node->getAttribute('op');
      my $att_providerId=$node->getAttribute('providerId');
      my $att_id_string=undef;
      my $att_op_string=undef;
      my $att_providerId_string = "";

      if (defined $att_id && $att_id ne ''){
         $att_id_string=" id=\"". $att_id."\"";
      }
      if (defined $att_op && $att_op ne ''){
         $att_op_string=" op=\"". $att_op."\"";
      }
      if (defined $att_providerId && $att_providerId ne ''){
         $att_providerId_string=" providerId=\"". $att_providerId."\"";
      }

      #warn "\nnode:", $node->getNodeName();
      print OUT "\n", " "x $indent, "<",  $node->getNodeName,  $att_providerId_string, ">";
      foreach my $child ($node->getChildNodes()) {
        _traverse($child, $indent+3);
      }
      if ($node->getFirstChild()->getNodeType() ==TEXT_NODE){
          print OUT "</", $node->getNodeName, ">";
      }
      else{
         print OUT "\n", " "x $indent,"</", $node->getNodeName, ">";
      }
    } elsif ($node->getNodeType() == TEXT_NODE) {
          my $str=$node->getData;
      #    $str=&_decon_xml($str);
          # print $node->getData;
          print OUT $str;
    }
  }

