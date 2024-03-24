#!/usr/local/bin/perl -w 

=head1 NAME

 Utils.pm - a grab bag of useful functions

=head1 DESCRIPTION

Includes functions to retrieve commonly needed bits of info from the database
as well as do some useful sorts and other types of things

=head2 Methods

=over 12

=item C<get_organism_id>

 arg1 - database handle
 arg2 - string 'genus species'

 @return int organism_id

 NOTE: splits genus species on first space so can deal with things like
       'Drosophila pseudoobscura pseudoobscura'

=item C<get_organism_id_by_abbr>

 arg1 - database handle
 arg2 - string organism abbreviation

 @return int organism_id

=item C<get_genus_species>

 arg1 - database handle
 arg2 - int organism_id

 @return array array[0] = genus array[1] = species

=item C<get_genus_species_by_abbr>

 arg1 - database handle
 arg2 - string species abbreviation

 @return array array[0] = genus array[1] = species

=item C<get_organism_abbr>

 arg1 - database handle
 arg2 - int organism_id

 @return string organism abbreviation

=item C<get_type_id>

 arg1 - database handle
 arg2 - string cvterm name for feature type

 @return int type_id

 NOTE: only looks for terms in SO cv

=item C<get_type_name>

 arg1 - database handle
 arg2 - int cvterm_id

 @return string cvterm.name

 NOTE: will work for any cvterm NOT specific to SO

=item C<get_feature_id>

 arg1 - database handle
 arg2 - string feature.uniquename
 arg3 - int feature.type_id
 arg4 - int feature.organism_id

=item C<retrieve_feature_id>

semi generic function to return a feature_id for a feature given various acceptable formats of identifying info

@params
db => database handle
feature => can be feature_id, uniquename, hashref to hash of feature info that includes either feature_id or uniquename key

@return int feature_id

=item C<guess_feature_id>

 arg1 - database handle
 arg2 - string feature.uniquename or feature_id or hashref to feature info

 @return int feature_id

=item C<get_feature_id_by_uname>

 arg1 - database handle
 arg2 - string feature.uniquename

 @return int feature_id

 NOTE: will return only one result but warns if more than
       one exist

=item C<get_expression_id>

arg1 - database handle
arg2 - string FBex expression uniquename

@return int expression_id

=item C<get_expression_id_by_md5>

arg1 - database handle
arg2 - string expression md5checksum

@return int expression_id

=item C<get_expression_md5>

arg1 - database handle
arg2 - expression_id

@return str MD5checksum

=item C<get_annotation_id_by_uname>

 arg1 - database handle
 arg2 - string feature.uniquename

 @return string current annotation ID

=item C<get_cv_id>

 arg1 - database handle
 arg2 - string cv.name

 @return int cv_id

=item C<get_cvterm_id>

 arg1 - database handle
 arg2 - string cvterm.name
 arg3 - int optional cv_id
 arg4 - int optional boolean for is_obsolete

 @return int cvterm_id

 NOTE: you don't need to specify a cv_id but if you don't
 and there are more than one term with the same name you will
 get only 1 id which could be the WRONG ONE

 NOTE: by default function ignores obsolete terms but to remove this constraint
       specify any positive integer as 4th parameter

=item C<get_cvterm_id_by_name_cv>

 arg1 - database handle
 arg2 - string cvterm.name
 arg3 - string cv.name
 arg4 - int optional boolean for is_obsolete

 @return int cvterm_id

=item C<get_cvterm_info

 arg1 - database handle
 arg2 - cvterm_id

 @return hashref to hash with all unique key cvterm info
         cvterm_id, name, is_obsolete, dbxref_id, accession, db, cv

=item C<get_db_id>

 arg1 - database handle
 arg2 - string db.name

 @return int db_id

=item C<get_db_name>

 arg1 - database handle
 arg2 - int db.db_id

 @return string db.name

=item C<get_accessions>

 returns accessions and versions for dbxrefs of a given feature_id
can specify a db_id and whether you want only current accession

 params
 db => database handle
 feature =>feature_id or uniquename
 db_id => optional db_id
 is_current => optional any true value
 
 @returns an arrayref to an array whose members are hashref - accession => acc, version => version

=item C<get_pub_info>

 arg1 - database handle
 arg2 - int pub_id

 @return string pub.uniquename

=item C<get_pub_id>

 arg1 - database handle
 arg2 - string pub.uniquename

 @return int pub_id

=item C<get_analysis_id>

 arg1 - database handle
 arg2 - string analysis.program
 arg3 - string analysis.sourcename
 arg4 - optional string analysis.programversion

 @return int analysis_id

=item C<get_analysis>

 arg1 - database handle
 arg2 - int analysis_id

 @return array (program, sourcename, programversion)

=item C<get_residues>

 arg1 - database handle
 arg2 - int feature_id

 @return string feature.residues

=item C<get_subsequence>

 arg1 - database handle
 arg2 - int feature_id (feature must have residues)
 arg3 - fmin
 arg4 - fmax
 arg5 - optional boolean if true returns the reverse complement

 @return string subsequence of feature

=item C<get_featureprop_id>

 **NOTE: this should really be named get_featureprop_type_id

 arg1 - database handle
 arg2 - string name of type of featureprop
 arg3 - string optional cv.name default = 'annotation property type'

 @return int featureprop.type_id

=item C<get_feature_info>

 returns all unique key info for a feature in form that can be used to create chxml element
 given one or more pieces that will uniquely identify the feature
 if the elements provided do not provide enough unique info then warning and return undef

 params
 db => database handle
 feature_id => db feature_id
 uniquename => uniquename for feature
 organism_id => 
 organism => string for organism eg. 'Drosophila melanogaster'
 genus => genus must be accompanied by species
 species => must accompany genus
 type_id => feature.type_id
 type => name of feature type must be from SO
 
 @returns a hash with the following info
 feature_id => db id of the feature
 uniquename => for feature
 organism_id => for feature
 genus => 
 species =>
 type_id =>
 type => name of the feature type

=item C<get_feature_info_by_name_type>

returns all unique key info for a feature in form that can be used to create chxml element
 given a feature.name and feature.type name

arg1 - object dbh
arg2 - string feature.name
arg3 - string feature.type name

@returns a hash with the following info
 feature_id => db id of the feature
 uniquename => for feature
 organism_id => for feature
 genus => 
 species =>
 type_id =>
 type => name of the feature type

=item C<get_featureloc>

 returns featureloc info for a feature in form that can be used to create chxml element
 
 params
 db => database handle
 feature => feature_id or uniquename
 src => srcfeature_id or src uniquename (chromosome) - optional
 srctype => either chromosome_arm or golden_path_region or cvterm_id for these - optional
 rank => optional integer
  
 @returns a hashref with the following info
 srcfeature_id => db id for src
 strand => for feature
 fmin => relative to srcfeature_id
 fmax => relative to srcfeature_id

=item C<get_all_features_of_type>

 NOTE: will return info on ALL features of the specified type

 arg1 - database handle
 arg2 - string feature type name or int feature.type_id

 @return hashref to hash with feature_ids as keys and ref to hash with uniquekey info 
         about feature (uniquename, organism.genus, organism.species, type name) as values

=item C<get_subject_ids>

 returns an array of feature_ids for all features with child relationships
 to the provided feature_id, can be restricted based on relationship type
 and child feature type

 params
 db => database handle
 object => parent feature_id, uniquename or feature object hashref with ukey info
 reltype => relationship type optional
 ftype => feature type for child

 @return array of child (subject) feature_ids

=item C<get_object_ids>

 returns an array of feature_ids for all features with parent relationships
 to the provided feature_id, can be restricted based on relationship type
 and parent feature type

 params
 db => database handle
 subject => child feature_id, uniquename or feature object hashref with ukey info
 reltype => relationship type optional
 ftype => relationship type for parent

 @return array of child (subject) feature_ids


=item C<get_scaffolds>

 returns hashref that including info including residues of all scaffolds 
 of the specified type from the optional specified organism

=item C<get_etype_names>

 NOTE: fairly AEVAST specific

 arg1 - database handle
 arg2 - optional list of type_ids to look for

 @return hashref to hash with type_ids as keys and names as values

=item C<get_max_annotation_id>

 arg1 - database handle
 arg2 - organism id

 @return the maximum gene annotation id for that organism

=item C<is_earlier_than>

 arg1 - string date in format yyyy-mm-dd and optional hour:min:sec
 arg2 - string date in format yyyy-mm-dd and optional hour:min:sec

 NOTE: assumes 12:00:00.000001 if no timestamp provided

 @return 1 if arg1 is earlier than arg2 otherwise returns 0

=item C<trim>

 arg - string or array of strings
 @return string or array with all leading and trailing whitespace removed from string(s)

=item C<_all_true>

 arg - array
 @return true if all elements in the array are defined

=item C<get_date>

 @returns current date in yyyy-mm-dd format

=item C<create_id_string>

 NOTE: not too generic USE WITH CAUTION
 creates a string of feature_ids to use in select statements to limit
 queries to only specified features (but only for feature table aliased as f3)

 arg1 - array ref to a list of feature_ids
 @return - string formatted as 'and f3.feature_id IN (...);

=item C<_create_hashref>

 arg1 - arrayref of keys
 arg2..n = scalar values
 @return reference to created hash

=item C<arrays_r_identical>

 arg1 - arrayref to array1
 arg2 - arrayref to array2

 @return 1 if arrays are identical

=item C<by_srcfeature_id_fmin>

=item C<by_fmin>

=item C<by_arm_strand_min>

 sorts by location by sourcename strand and fmin

=item C<by_arm_min>

 sort by arm and min irrespective of strand

=item C<by_min>

 sort by min

=item C<by_max>

 sort by max

=item C<by_5to3>

 sort by 5' to 3' regardless of strand

=item C<reverse_complement>

 arg - string of nucleotide residues or array of strings
 @return reverse complemented string or array of strings

=item C<is_protein>

 arg - string of residues
 @return 1 if all chars are aminoacid symbols

=item C<is_protein_ambig>

 arg - string of residues
 @return 1 if all chars are within IUPAC ambiguity symbol table

=item C<is_rna>

 arg - string of residues
 @return 1 if residues are consistent with RNA i.e. no T's

=item C<is_dna>

 arg - string of residues
 @return 1 if residues are consistent with DNA i.e. ATGC

=item C<is_dna_ambig>

 arg - string of residues
 @return 1 if residues are consistent with DNA including N's

=item C<is_nucleic>

 arg - string of residues
 @return 1 if residues are consistent with any nucleic acid including ambiguity codes

=item C<get_formatted_time>

 Gets localtime() and returns the time stamp in this format: YYYY-MM-DD HH:MM:SS.

=item C<get_highest_caller_line_number>

 Gets the highest stack level from which a function was called and returns the line number.
 This is useful in logging.

=item C<print_log>

 arg - text to print out
 Adds timestamp to start of print text; makes debug of slow steps easier.

=head1 AUTHOR

Andy Schroeder - andy@morgan.harvard.edu
Gil dos Santos - dossantos@morgan.harvard.edu

=cut


#############################################################################################
######################### general utility functions #########################################
#############################################################################################
# return the organism_id for the given feature
sub get_organism_id_for_feature {
  my $dbh = shift;
  my $feature = shift;
  my $fid = guess_feature_id($dbh, $feature);
  return $dbh->selectrow_array("SELECT organism_id FROM feature WHERE feature_id = $fid");
}

# will determine the organism_id for the critter specified by genus species string
sub get_organism_id {
    my $dbh = shift;
    my $organism = shift;
    my ($genus, $species) = split (/\s+/,$organism,2);
    $species = '' unless $species;
    my ($org_id) = $dbh->selectrow_array
      (sprintf("SELECT organism_id FROM organism WHERE genus = '$genus' and species = '$species'"));
    print "WARNING: Can't find organism_id for $organism\n" and return unless $org_id;
    return $org_id;
}

# will determine the organism_id for the critter specified by abbreviation
sub get_organism_id_by_abbr {
  my $dbh = shift;
  my $abbr = shift;
  (my $cnt) = $dbh->selectrow_array(sprintf("SELECT count(*) FROM organism WHERE abbreviation = '$abbr'"));
  if ($cnt < 1) {
    print "WARNING: Can't find organism_id for $abbr\n" and return;
  } else {
    print "WARNING: More than one organism are abbreviated $abbr\n\tONLY THE FIRST RETRIEVED IS RETURNED"
      if $cnt > 1;
    (my $oid) = $dbh->selectrow_array(sprintf("SELECT organism_id FROM organism WHERE abbreviation = '$abbr'"));
    return $oid;
  }
}
  
# will return genus and species for specified organism_id
sub get_genus_species {
    my $dbh = shift;
    my $org_id = shift;
    my @genspec = $dbh->selectrow_array
      (sprintf("SELECT genus, species FROM organism WHERE organism_id = $org_id"));
    print "WARNING: Can't find genus species info for organism_id = $org_id\n" and return unless @genspec;
    return @genspec;
}

# returns genus and species in that order for given abbreviation
sub get_genus_species_by_abbr {
  my $dbh = shift;
  my $abbr = shift;
  my $oid = get_organism_id_by_abbr($dbh,$abbr);
  return unless $oid;
  return get_genus_species($dbh,$oid);
}

sub get_organism_abbr {
  my $dbh = shift;
  my $oid = shift;
  my ($abbr) = $dbh->selectrow_array
      (sprintf("SELECT abbreviation FROM organism WHERE organism_id = $oid"));
  print "WARNING: Can't find abbreviation for organism with id $oid\n" and return unless $abbr;
  return $abbr;
}

sub get_all_so_terms {
  my $dbh = shift;
  my %types;
  my $stmt = "SELECT cvterm_id, c.name FROM cvterm c, cv WHERE c.cv_id = cv.cv_id and cv.name = 'SO' and is_obsolete = 0";
  my $q = $dbh->prepare($stmt); 
  $q->execute or die "Can't get SO terms for feature types\n";
  while ((my $id, my $term) = $q->fetchrow_array()) {
    $types{$term} = $id;
  }
  return %types;
}


#specific for SO feature types
sub get_type_id {
     my $dbh = shift;
     my $term = shift;
     my ($type_id) = $dbh->selectrow_array
       ("SELECT cvterm_id FROM cvterm c, cv
         WHERE c.cv_id = cv.cv_id and cv.name = 'SO' and c.name = '$term' and c.is_obsolete = 0");
     print "NO TYPE ID WAS FOUND FOR $term\n" and return unless $type_id;
     return $type_id;
}

# returns the cvterm.name for the specified type_id/cvterm_id
sub get_type_name {
     my $dbh = shift;
     my $type_id = shift;
     my ($name) = $dbh->selectrow_array
       ("SELECT name FROM cvterm WHERE cvterm_id = $type_id");
     print "CAN'T FIND A NAME FOR TYPE_ID = $type" and return unless $name;
     return $name;
}

sub guess_feature_id {
  my $dbh = shift;
  my $feat = shift;
  my $feature_id;
  if ($feat =~ /^[0-9]+$/) {
    # assume a feature_id
    $feature_id = $feat;
  } elsif (ref($feat) eq 'HASH') {
    if ($feat->{feature_id}) {
      $feature_id = $feat->{feature_id};
    } else {
      $feat->{db} = $dbh;
      my $nfeat = get_feature_info($feat);
      $feature_id = $nfeat->{feature_id};
    }
  } else {
    # assume a uniquename
    $feature_id = get_feature_id_by_uname($dbh, $feat);
  }

  print STDERR "CAN'T GUESS THE feature_id -- NO GO!\n" and return unless $feature_id;

  if (@_) {
    unless ($dbh->selectrow_array("SELECT feature_id FROM feature WHERE feature_id = $feature_id")) {
      print STDERR "WARNING -- CAN'T VERIFY THE FEATURE ID AS VALID - RETURNING!\n"
	and return;
    }
  }
  return $feature_id;
}

# given uniquename, type_id and organism_id will return feature_id
sub get_feature_id {
    my $dbh = shift;
    my $uname = shift;
    my $tid = shift;
    my $oid = shift;
    my $is_curr = shift if @_;

    my $stmt = "SELECT feature_id FROM feature WHERE uniquename = '$uname'
        and    type_id = $tid and organism_id = $oid";
    $stmt .= " and is_obsolete = false" if ($is_curr); 
    my ($fid) = $dbh->selectrow_array($stmt);
    print "Can't find the feature with uniquename $uname\n" and return unless $fid;
    return $fid;
}

sub get_feature_id_by_uname {
  my $dbh = shift;
  my $uname = shift;

  my ($cnt) = $dbh->selectrow_array
    (sprintf
     ("SELECT count(*) FROM feature WHERE uniquename = '$uname' and is_obsolete = false"));
  if ($cnt < 1) {
    print "WARNING: NO FEATURE WITH $uname FOUND\n" and return;
  } else {
    print "WARNING: More than one feature with $uname FOUND!\n\tONLY THE FIRST RETRIEVED IS RETURNED" if $cnt > 1;
    return $dbh->selectrow_array(sprintf("SELECT feature_id FROM feature WHERE uniquename = '$uname' and is_obsolete = false"));
  }
}

# given expression uniquename returns expression_id
sub get_expression_id {
  my $dbh = shift;
  my $uname = shift;

  return $dbh->selectrow_array(sprintf("SELECT expression_id FROM expression WHERE uniquename = '$uname'"));
}

# given expression md5checksum returns expression_id
sub get_expression_id_by_md5 {
  my $dbh = shift;
  my $md5 = shift;

  return $dbh->selectrow_array(sprintf("SELECT expression_id FROM expression WHERE md5checksum = '$md5'"));
}

# given expression_id returns md5checksum
sub get_expression_md5 {
  my $dbh = shift;
  my $eid = shift;

  return $dbh->selectrow_array(sprintf("SELECT md5checksum FROM expression WHERE expression_id = $eid"));
}



# given uniquename returns the current annotation_id for the feature
sub get_annotation_id_by_uname {
  my $dbh = shift;
  my $uname = shift;

  (my $cnt) = $dbh->selectrow_array(sprintf("SELECT count(*) FROM dbxref d, feature_dbxref fd, feature f, db
                                             WHERE  d.dbxref_id = fd.dbxref_id and d.db_id = db.db_id
                                             and    fd.feature_id = f.feature_id and fd.is_current = true
                                             and    db.name = 'FlyBase Annotation IDs'
                                             and    f.uniquename = '$uname'"));
  if ($cnt < 1) {
    print "WARNING: NO ANNOTATION ID FOUND FOR $uname\n" and return;
  } else {
    print "WARNING: More than one annotation id was found for $uname !\n\tONLY THE FIRST RETRIEVED IS RETURNED"
      if $cnt > 1;
    (my $aid) = $dbh->selectrow_array(sprintf("SELECT accession FROM dbxref d, feature_dbxref fd, feature f, db
                                             WHERE  d.dbxref_id = fd.dbxref_id and d.db_id = db.db_id
                                             and    fd.feature_id = f.feature_id and fd.is_current = true
                                             and    db.name = 'FlyBase Annotation IDs'
                                             and    f.uniquename = '$uname'"));
    return $aid;
  }
}



sub get_cv_id {
    my $dbh = shift;
    my $cv = shift;

    my ($cvid) = $dbh->selectrow_array
      ("SELECT cv_id FROM cv WHERE name = '$cv'");
    print "Can't find cv name $cv\n" and return unless $cvid;
    return $cvid;
}

# will return a cvterm_id for a cvterm
# WARNING you don't need to specify a cv_id but if you don't
# and there are more than one term with the same name you will
# get only 1 id which could be the WRONG ONE
# by default ignores obsolete terms but to remove this constraint
# specify an integer as 4th parameter
sub get_cvterm_id {
    my $dbh = shift;
    my $name = shift;
    my $cv_id;
    $cv_id = shift if @_;
    my $cv_id_string = '';
    $cv_id_string = "and cv_id = $cv_id" if $cv_id;
    $obs = shift if @_;
    my $obs_string = 'and is_obsolete = 0';
    $obs_string = '' if $obs;

    my ($cvtid) = $dbh->selectrow_array
      ("SELECT cvterm_id FROM cvterm WHERE name = '$name' $cv_id_string $obs_string");
    print "Can't find cvterm name $name\n" and return unless $cvtid;
    return $cvtid;
}

sub get_cvterm_id_by_name_cv {
  my $dbh = shift;
  my $name = shift;
  my $cv = shift;
  my $is_obs = shift if @_;

  my $obs_string = ' and is_obsolete = 0';
  $obs_string = '' if $is_obs;

  (my $cvtid) = $dbh->selectrow_array
    (sprintf
     ("SELECT cvterm_id FROM cvterm c, cv
      WHERE  c.name = '$name' and c.cv_id = cv.cv_id
        and  cv.name = '$cv' $obs_string"));
  print "Can't find cvterm name $name in CV $cv\n" 
    and return unless $cvtid;
  return $cvtid;
}

# given a cvterm_id get all unique key info
sub get_cvterm_info {
  my $dbh = shift;
  my $cvtid = shift;

  my @cvterm = $dbh->selectrow_array
    (sprintf
     ("SELECT c.name, c.is_obsolete, d.dbxref_id, d.accession, db.name, cv.name, c.cvterm_id
       FROM   cvterm c, dbxref d, cv, db
       WHERE  c.cvterm_id = $cvtid
         and  c.dbxref_id = d.dbxref_id and d.db_id = db.db_id
         and  c.cv_id = cv.cv_id"
     )
    );
  return {name => $cvterm[0],
	  is_obsolete => $cvterm[1],
	  dbxref_id => $cvterm[2],
	  accession => $cvterm[3],
	  db => $cvterm[4],
	  cv => $cvterm[5],
          cvterm_id => $cvterm[6],
	 };
}

# specify a db name and it will return the db_id
sub get_db_id {
  my $dbh = shift;
  my $name = shift;
  (my $db_id) = $dbh->selectrow_array("SELECT db_id from db where name = '$name'");
  print "NO DB with name:$name exists in the instance you are using\n"
    and return unless $db_id;
  return $db_id;
}

# specify a db_id and it will return the db.name
sub get_db_name {
  my $dbh = shift;
  my $id = shift;
  (my $name) = $dbh->selectrow_array("SELECT name from db where db_id = $id");
  print "NO DB with name:$name exists in the instance you are using\n"
    and return unless $name;
  return $name;
}

# return accessions (and versions) given a feature_id and one or more optional params
sub get_accessions {
  my %params = @_;
  print "ERROR - missing database handle!\n" and return unless $params{db};
  print "ERROR - missing feature id!\n" and return unless $params{feature};
  my $dbh = $params{db};
  my $fid = guess_feature_id($dbh, $params{feature});
  my $dbid;
  my $is_current;
  $dbid = $params{db_id} if $params{db_id};
  $is_current = 'true' if $params{is_current};

  my $stmt = "SELECT DISTINCT accession, version FROM feature_dbxref fd, dbxref d WHERE fd.dbxref_id = d.dbxref_id and fd.feature_id = $fid";

  $stmt .= " and d.db_id = $dbid" if $dbid;
  $stmt .= " and fd.is_current = true" if $is_current;

  my $acc_q = $dbh->prepare($stmt);
  $acc_q->execute or die "Can't query for accessions\n";

  my @accs;
  while ((my $acc, $version) = $acc_q->fetchrow_array()) {
    push @accs, {accession => $acc, version => $version};
  }
  return \@accs if @accs;
  return;
}


# returns the pub uniquename (at the moment the unique key on pub)
# if provided a valid pub_id
sub get_pub_info {
    my $dbh = shift;
    my $pid = shift;

    my ($puname) = $dbh->selectrow_array
      ("SELECT uniquename FROM pub WHERE pub_id =  $pid");
    print "Can't find pub uniquename $puname\n" and return unless $puname;
    return $puname;
}

# returns the pub id (at the moment the unique key on pub)
# if provided a valid pub uniquename
sub get_pub_id {
    my $dbh = shift;
    my $puname = shift;

    my ($pid) = $dbh->selectrow_array
      ("SELECT pub_id FROM pub WHERE uniquename = '$puname'");
    print "Can't find pub id for $puname\n" and return unless $pid;
    return $pid;
}

# given a feature_id will return the residues of that feature
# arg1 = database handle, arg2 = feature_id
sub get_residues {
    my $dbh = shift;
    my $fid = shift;
    my ($res) = $dbh->selectrow_array
      ("SELECT residues FROM feature WHERE feature_id = $fid");
    return $res;
}

sub get_subsequence {
  my $dbh = shift;
  my $fid = shift;
  my $fmin = shift;
  my $fmax = shift;
  my $reverse = shift;

  (my $seq) = $dbh->selectrow_array
      ("SELECT substring(residues FROM $fmin + 1 FOR $fmax - $fmin) FROM feature WHERE feature_id = $fid");
  return unless $seq;
  return reverse_complement($seq) if $reverse;
  return $seq;
}

# given program and sourcename (and optional program version) returns analysis id
sub get_analysis_id {
  my $dbh = shift;
  my $prog = shift;
  my $srcname = shift;
  my $vers = shift if @_;

  my $stmt = "SELECT analysis_id FROM analysis WHERE program = '$prog' and sourcename = '$srcname'";
  $stmt .= " and programversion = $vers" if $vers;

  (my $analysis_id) = $dbh->selectrow_array($stmt);
  return $analysis_id;
}

# given analysis_id return program, sourcename and programversion
sub get_analysis {
  my $dbh = shift;
  my $aid = shift;

  my $stmt = "SELECT program, sourcename, programversion FROM analysis WHERE analysis_id = $aid";
  (my @analysis) = $dbh->selectrow_array($stmt);
  return @analysis;
}
  

# returns the featureprop_id for the specified featureprop from cv specified
# note that if no cv is specified will default to 'annotation property type' cv
# will only return a single id
# 
# arg1 = db handle
# arg2 = name of term you want to get id for
# optional arg3 = name of cv to get the term from
sub get_featureprop_id {
     my $dbh = shift;
     my $term = shift;
     my $cv = 'annotation property type'; #default cv
     $cv = shift if @_; #change cv if one is provided
     my ($type_id) = $dbh->selectrow_array
       ("SELECT cvterm_id FROM cvterm c, cv
         WHERE c.cv_id = cv.cv_id and cv.name = '$cv' and c.name = '$term'");
     print "NO TYPE ID WAS FOUND FOR $term from CV $cv\n" and return unless $type_id;
     return $type_id;
}

# returns all unique key info for a feature in form that can be used to create chxml element
# given one or more pieces that will uniquely identify the feature
# if the elements provided do not provide enough unique info then warning and return undef
#
# params
# db => database handle
# feature_id => db feature_id
# uniquename => uniquename for feature
# organism_id => 
# organism => string for organism eg. 'Drosophila melanogaster'
# genus => genus must be accompanied by species
# species => must accompany genus
# type_id => feature.type_id
# type => name of feature type must be from SO
# analysis => true or false if you want to include this filter
# 
# returns a hash with the following info
# feature_id => db id of the feature
# uniquename => for feature
# organism_id => for feature
# genus => 
# species =>
# type_id =>
# type => name of the feature type
sub get_feature_info {
  my %params = @_;
  print "WARNING - no database handle provided - NO GO!\n" and return unless $params{db};
  my $dbh = $params{db};

  # here is the base statement - add various bits depending what info is passed
  my $stmt = "SELECT feature_id, uniquename, f.organism_id, o.genus, o.species, f.type_id, c.name FROM feature f, organism o, cvterm c WHERE f.organism_id = o.organism_id and f.type_id = c.cvterm_id and f.is_obsolete = false";

  if ($params{analysis}) {
    my $analysis;
    if ($params{analysis} == 1 or $params{analysis} =~ /t/i) {
      $analysis = 'true';
    } else {
      $analysis = 'false';
    }
    $stmt .= " and is_analysis = $analysis";
  }

  if ($params{feature_id}) { # this is all the info we need
    $stmt .= " and feature_id = $params{feature_id}";
  } else {
    $stmt .= " and uniquename = '$params{uniquename}'" if $params{uniquename};

    # get organism_id if any organism params passed
    my $org_id;
    if ($params{organism_id}) {
      $org_id = $params{organism_id};
    } elsif ($params{organism}) {
      $org_id = get_organism_id($dbh, $params{organism});
    } elsif ($params{genus} and $params{species}) {
      my $orgn = "$params{genus} $params{species}";
      $org_id = get_organism_id($dbh, $orgn);
    }
    $stmt .= " and o.organism_id = $org_id" if $org_id;

    # deal with type
    my $tid;
    if ($params{type_id}) {
      $tid = $params{type_id};
    } elsif ($params{type}) {
      $tid = get_type_id($dbh,$params{type});
    }
    $stmt .= " and type_id = $tid" if $tid;
  }

  #print "PREPARING TO EXECUTE $stmt\n";
  my $cntstmt = $stmt;
  $cntstmt =~ s/SELECT.*FROM/SELECT COUNT(1) FROM/;
  (my $cnt) = $dbh->selectrow_array($cntstmt);
  if ($cnt < 1) {
    print "WARNING: can't find a feature based on the provided info - NO GO!\n" and return;
  } elsif ($cnt > 1) {
    print "WARNING: more than one feature are identified based on the provided info - NO GO!\n" and return;
  } else {
    my @info = $dbh->selectrow_array($stmt);
    return (feature_id => $info[0],
	    uniquename => $info[1],
	    organism_id => $info[2],
	    genus => $info[3],
	    species => $info[4],
	    type_id => $info[5],
	    type => $info[6],
	   );
  }
}

# provide a feature.name and feature.type (name) return info about the feature if found
sub get_feature_info_by_name_type {
  my $dbh = shift;
  my $feat = shift;
  my $type = shift;
  my $alive = shift;

  #print "RETRIEVING INFO FOR $feat of TYPE $type\n";

  my $tyid = get_type_id($dbh, $type);
  unless ($tyid) {
    print "WARNING - couldn't find $type in the DB\n";
    return;
  }
  my $stmt = "SELECT feature_id FROM feature WHERE name = '$feat' and type_id = $tyid";

  $stmt .= " and is_obsolete = false" if $alive;

  my $f_q = $dbh->prepare($stmt);
  $f_q->execute or die "Can't execute $stmt\n";
  if ($f_q->rows == 0) {
    warn "WARNING - NO feature of name $feat and type $type FOUND - returning!\n";
    return;
  } elsif ($f_q->rows > 1) {
    unless ($alive) {
      return get_feature_info_by_name_type($dbh, $feat, $type, 1);
    } else {
      warn  "WARNING - MORE THAN ONE feature of name $feat and type $type FOUND - returning!\n";
      return;
    }
  } 
  (my $fid) = $f_q->fetchrow_array();
  #print "CHECKING FEATURE_ID $fid\n";

  my %finfo = get_feature_info(db => $dbh, feature_id => $fid);
  return \%finfo;
}

  
# returns a hashref to featureloc info
# fmin, fmax, srcfeature_id, strand for a feature
sub get_featureloc {
  my %params = @_;
  print "WARNING -- missing required parameter database handle or feature info - NO GO!\n"
    and return unless ($params{db} and $params{feature});

  my $dbh = $params{db};
  my $feat = $params{feature};
  my $fid = guess_feature_id($dbh, $feat);
  print "WARNING - Couldn't get a feature_id -- RETURNING\n"
    and return unless $fid;

  # look for optional parameters
  my $srcstr = ''; my $rankstr = '';

  if ($params{src}) {
    my $srcid = guess_feature_id($dbh, $params{src});
    $srcstr = " and srcfeature_id = $srcid " if $srcid;
  }

  my $srctyid;
  if (my $srctype =$params{srctype}) {
    if ($srctype =~ /^[0-9]+$/) {
      $srctyid = $srctype;
    } else {
      $srctyid = get_type_id($dbh, $srctype);
    }
  }
  my $ftablestr = '';
  my $srctypestr = '';
  if ($srctyid) {
    $ftablestr = ' , feature ';
    $srctypestr = " and srcfeature_id = feature.feature_id and feature.type_id = $srctyid ";
  }

  if ($params{rank}) {
    $rankstr = " and rank = $params{rank} " if ($params{rank} =~ /^[0-9]$/);
  }

  my $stmt = sprintf
     ("SELECT featureloc.*
       FROM   featureloc %s
       WHERE  featureloc.feature_id = $fid %s %s %s", 
      $ftablestr, $srcstr, $rankstr, $srctypestr
     );
#  print "$stmt\n";
  my $fl_q = $dbh->prepare($stmt);

  $fl_q->execute or die "Can't do location query\n";

  print "WARNING - more than one row retrieved - returning first row only\n" 
    if $fl_q->rows > 1;

  return $fl_q->fetchrow_hashref;
}


# arg1 = database handle, arg2 = feature.type name or id
# will return a hashref of all features of the specified type
# hash will have feature_id as key and ref to hash with uniquekey info 
# about feature (uniquename, organism.genus, organism.species, type name) as value
#
# arg1 = database handle, arg2 = feature.type_id
sub get_all_features_of_type {
    my $dbh = shift;
    my $tid = shift;
    $tid = get_type_id($dbh,$tid) if ($tid !~ /^[0-9]+$/); # a type name has been specified so get type_id

    my %features;

    print "Getting all features of type $tid\n";
    my $query = $dbh->prepare
      (sprintf
       ("SELECT feature_id, uniquename, genus, species, c.name 
         FROM   feature f
         JOIN   organism o ON (f.organism_id = o.organism_id and f.type_id = $tid)
         JOIN   cvterm c ON (f.type_id = c.cvterm_id) and is_obsolete = false"
       )
      );

    $query->execute or die "Can't retrieve features of type id $tid\n";

    while (my ($fid, $uname, $genus, $species, $type) = $query->fetchrow_array()) {
      $features{$fid} = {
			 uniquename => $uname,
			 genus => $genus,
			 species => $species,
			 type => $type,
			};
    }

    %features;
}

# when provided a feature_id and optional feature_relationship type and 
# child feature type returns an array of subject features
sub get_subject_ids {
  my %params = @_;
  print "No database handle - NO GO!\n" and return unless $params{db};
  print "No object feature_id provided - NO GO!\n" and return unless $params{object};

  my $dbh = $params{db};
  my $obj = $params{object};
  # figure out which of the options was provided as obj and translate into a feature_id
  my $object_id = guess_feature_id($dbh,$obj);

  print "Can't find a valid id for the object you provided - NO GO!\n"
    and return unless $object_id;

  # check for optional params
  my $relstring = '';
  if ($params{reltype}) {
    my $reltype = $params{reltype};
    my $rtype_id;
    if ($reltype =~ /^[0-9]+$/) {
      $rtype_id = $reltype;
    } else {
      my $cv_id = get_cv_id($dbh, 'relationship type');
      $rtype_id = get_cvterm_id($dbh, $reltype, $cv_id);
    }
    if ($rtype_id) {
      $relstring = " and fr.type_id = $rtype_id";
    } else {
      print "Can't find a relationship of type $reltype - will return all relationship types\n";
    }
  }

  my $ftable_string = '';
  my $ftype_string = '';

  if ($params{ftype}) {
    my $tid;
    my $type = $params{ftype};
    if ($type =~ /^[0-9]+$/) {
      # assume valid feature type
      $tid = $type;
    } else {
      # assume type name
      my $cv_id = get_cv_id($dbh, 'SO');
      $tid = get_cvterm_id($dbh, $type, $cv_id);
    }
    if ($tid) {
      $ftable_string = ", feature f ";
      $ftype_string = " and f.feature_id = fr.subject_id and f.type_id = $tid ";
    } else {
      print "Can't find a valid type_id for $type so fetching all related feature types\n";
    }
  }

  my $query = $dbh->prepare
    (sprintf
     ("SELECT fr.subject_id FROM feature_relationship fr %s
       WHERE  fr.object_id = $object_id %s %s", $ftable_string, $relstring, $ftype_string)
    );

  $query->execute or die "Can't execute child query\n";

  my @subjects;
  while ((my $fid) = $query->fetchrow_array) {
    push @subjects, $fid;
  }
  return @subjects;
}

# when provided a feature_id and optional feature_relationship type and 
# child feature type returns an array of object features
sub get_object_ids {
  my %params = @_;
  print "No database handle - NO GO!\n" and return unless $params{db};
  print "No object feature_id provided - NO GO!\n" and return unless $params{subject};

  my $dbh = $params{db};
  my $obj = $params{subject};
  # figure out which of the options was provided as obj and translate into a feature_id
  my $object_id = guess_feature_id($dbh,$obj);

  print "Can't find a valid id for the subject you provided - NO GO!\n"
    and return unless $object_id;

  # check for optional params
  my $relstring = '';
  if ($params{reltype}) {
    my $reltype = $params{reltype};
    my $rtype_id;
    if ($reltype =~ /^[0-9]+$/) {
      $rtype_id = $reltype;
    } else {
      my $cv_id = get_cv_id($dbh, 'relationship type');
      $rtype_id = get_cvterm_id($dbh, $reltype, $cv_id);
    }
    if ($rtype_id) {
      $relstring = " and fr.type_id = $rtype_id";
    } else {
      print "Can't find a relationship of type $reltype - will return all relationship types\n";
    }
  }

  my $ftable_string = '';
  my $ftype_string = '';

  if ($params{ftype}) {
    my $tid;
    my $type = $params{ftype};
    if ($type =~ /^[0-9]+$/) {
      # assume valid feature type
      $tid = $type;
    } else {
      # assume type name
      my $cv_id = get_cv_id($dbh, 'SO');
      $tid = get_cvterm_id($dbh, $type, $cv_id);
    }
    if ($tid) {
      $ftable_string = ", feature f ";
      $ftype_string = " and f.feature_id = fr.object_id and f.type_id = $tid ";
    } else {
      print "Can't find a valid type_id for $type so fetching all related feature types\n";
    }
  }

  my $query = $dbh->prepare
    (sprintf
     ("SELECT fr.object_id FROM feature_relationship fr %s
       WHERE  fr.subject_id = $object_id %s %s", $ftable_string, $relstring, $ftype_string)
    );

  $query->execute or die "Can't execute child query\n";

  my @subjects;
  while ((my $fid) = $query->fetchrow_array) {
    push @subjects, $fid;
  }
  return @subjects;
}


sub get_scaffolds {
  my %params = @_;

  print "NO DATABASE HANDLE PROVIDED -- NO GO!\n" and return unless $params{db};
  my $dbh = $params{db};

  print "YOU MUST SPECIFY A FEATURE TYPE TO RETRIEVE (chromosome_arm, golden_path_region or a type_id) -- NO GO!\n" and return unless $params{type};

  my $type = $params{type};

  my $tid;

  if ($type =~ /^[0-9]+$/) {
    # assume a valid type_id is provided
    $tid = $type;
  } else {
    $tid = get_type_id($dbh,$type);
  }

  my $orgn_str = '';
  if ($params{organism}) {
    my $orgn = $params{organism};
    my $oid;
    if ($orgn =~ /^[0-9]+$/) {
      $oid = $orgn;
    } else {
      $oid = get_organism_id($dbh,$orgn);
    }
    $orgn_str = " and organism_id = $oid ";
  }

  my $sc_q = $dbh->prepare
    (sprintf
     ("SELECT feature_id, uniquename, residues, organism_id
       FROM   feature WHERE type_id = $tid and is_obsolete = false
       %s", $orgn_str));

  $sc_q->execute or die "Can't do scaffold query\n";

  my %scaffs;
  while ((my $fid, my $uname, my $res, my $oid) = $sc_q->fetchrow_array) {
    $scaffs{$fid} = {feature_id => $fid,
		     uniquename => $uname,
		     residues => $res,
		     organism_id => $oid,};
  }
  return \%scaffs;
}

sub get_feature_id_by_accession {
  my %params = @_;
  my $dbh = $params{dbh};
  my $acc = $params{acc};
  my $version;
  my $db;
  my $is_curr = 1;
  my $ftype;
  warn "Missing required info\n" and return unless $dbh and $acc;
  $version = $params{version} if $params{version};
  $db = $params{db} if $params{db};
  if ($params{curr}) {
    undef $is_curr unless ($params{curr} == 1 or $params{curr} =~ /^t/i);
  }
  $ftype = $params{ftype} if $params{ftype};
  my $ftyid = get_type_id($dbh, $ftype) if $ftype;

  my $tbl2add = ' ';
  my $conds2add = ' ';

  if ($db) {
    $tbl2add .= ', db ';
    $conds2add .= " and d.db_id = db.db_id and db.name = '$db' ";
  }

  $conds2add .= " and version = $version " if $version;
  $conds2add .= " and is_current = true " if $is_curr;
  $conds2add .= " and f.type_id = $ftyid " if $ftyid;

  my $stmt = sprintf( "SELECT f.feature_id FROM feature f, feature_dbxref fd, dbxref d %s
                       WHERE f.feature_id = fd.feature_id and fd.dbxref_id = d.dbxref_id and f.is_obsolete = false
                               and  accession = '$acc' %s", $tbl2add, $conds2add);

  my $query = $dbh->prepare($stmt);
  $query->execute or return "Can't execute $stmt\n";
  if ($query->rows == 0) {
    return;
  } elsif ($query->rows > 1) {
    print "WARNING more than one feature associated with $acc found ONLY returning the first\n";
  }
  (my $result) = $query->fetchrow_array();
  return $result;
}


# function to make a hash with evidence type_ids as key and names as hash
# optional arg - arrayref to array of types of evidence to look for overlap 
# return hashref to hash with type_ids as keys and names as values
sub get_etype_names {
    my $dbh = shift;
    my %type2text;
    my @type_ids;
    
    if (@_) {
	print "USING TYPE NAMES PROVIDED\n";
	@type_ids = @{+shift};
    } else {
	print "GETTING TYPE NAMES FROM DB\n";
	#query chado for all types of evidence
	my $query = $dbh->prepare
	    (sprintf
	     ("SELECT DISTINCT f2.type_id 
               FROM feature f1, featureloc fl, feature f2
               WHERE f1.is_analysis = 'true'
                 and f1.feature_id = fl.feature_id
                 and rank = 1
                 and srcfeature_id = f2.feature_id"));
	$query->execute or die "Cannot execute get_types query\n";
	while ((my $t_id) = $query->fetchrow_array()) {
	    push(@type_ids, $t_id);
	}
    }
    
    foreach my $type (@type_ids) {
	my $query = $dbh->prepare
	    (sprintf("SELECT name FROM cvterm WHERE cvterm_id = $type"));
	$query->execute or die "Cannot execute type2text query\n";
	
	while ((my $name) = $query->fetchrow_array()) {
	    $type2text{$type} = $name;
	}
    }
    return \%type2text;
}

# want to turn into a kind of generic fxn to return a feature_id
# given any variety of possible inputs
sub retrieve_feature_id {
  my %params = @_;
  print "Database handle and aligned feature required -- NO GO!\n"
    and return unless $params{db} and $params{feature};
  my $dbh = $params{db};
  my $feat = $params{feature};

  # for feature expect either a hashref with feature_id or uniquename
  # one of the values or a feature_id or uniquename and return otherwise
  my $fid;
  if (ref($feat) eq 'HASH') {
    if ($feat->{feature_id}) {
      $fid = $feat->{feature_id};
    } elsif ($feat->{uniquename}) {
      $fid = get_feature_id_by_uname($dbh, $feat->{uniquename});
    }
  } elsif ($feat =~ /^[0-9]+$/) {
    $fid = $feat;
  } else {
    $fid = get_feature_id_by_uname($dbh,$feat);
  }
  print "Couldn't identify a feature_id - NO GO!\n" and return unless $fid;
  return $fid;
}

# returns the max annotation id for the organism provided
sub get_max_annotation_id {
  my $dbh = shift;
  my $orgid = shift;

  unless ($orgid =~ /^[0-9]+$/) {
    if ($orgid =~ /^D[a-z]{3}$/) { # it's an abbreviation
      $orgid = get_organism_id_by_abbr($dbh, $orgid);
    } else { # assume genus species
      $orgid = get_organism_id($dbh, $orgid);
    }
  }

  print "Can't determine the organism\n" and return unless $orgid;

  (my $db_id) = $dbh->selectrow_array("SELECT db_id FROM db WHERE name = 'FlyBase Annotation IDs'");
  my $melid = get_organism_id_by_abbr($dbh,'Dmel');
  my $gtyid = get_type_id($dbh,'gene');

  return $dbh->selectrow_array("SELECT max(accession) FROM feature f, feature_dbxref fd, dbxref d WHERE f.feature_id = fd.feature_id and fd.dbxref_id = d.dbxref_id and f.organism_id = $orgid and f.type_id = $gtyid and d.db_id = $db_id and accession SIMILAR TO '(C|G)[A-Z][0-9]{5}'");

}

sub get_exons_of_transcript {
  my $dbh = shift;
  my $feat_id = shift;
  my $srcfeat_type = shift;
  my $exon_fl_rank = shift || 0;

  #print "srcfeat_type: $srcfeat_type\n";
  my $srcfeat_type_id = get_type_id($dbh, $srcfeat_type);
  my $exon_ftypeid = get_type_id($dbh, 'exon');

  #query for annotation gene models
  my $sql = "select l.fmin, l.fmax, l.strand, l.phase, f.uniquename from featureloc l, feature_relationship r, cvterm t, feature f, feature s where f.type_id = $exon_ftypeid and f.feature_id = l.feature_id and l.feature_id = r.subject_id and r.object_id = $feat_id and r.type_id = t.cvterm_id and t.name = 'partof' and l.srcfeature_id = s.feature_id and s.type_id = $srcfeat_type_id and l.rank = $exon_fl_rank order by r.rank";

  #query for gene predictions stored as evidence
  #my $sql = "select l.fmin, l.fmax, l.strand, l.phase, f.uniquename from featureloc l, feature_relationship r, cvterm t, feature f, feature s where f.feature_id = l.feature_id and l.feature_id = r.subject_id and r.object_id = $mrna_id and r.type_id = t.cvterm_id and t.name = 'partof' and l.srcfeature_id = s.feature_id and s.type_id = $srcfeat_type_id and l.rank = $exon_fl_rank order by l.fmin * l.strand";
  #print $sql, "\n";
  my $q = $dbh->prepare($sql);
  $q->execute || die $q->errstr;
  
  my $ref = $q->fetchall_arrayref;

  return $ref;
}

#the following sub retrieves the sequence of a transcript
sub get_transcript_seq {
  my $dbh = shift;
  my $tr_id = shift;
  my $srcfeat_type = shift;
  my $fl_rank = shift || 0;	
  my $is_prediction = shift or undef; # not currently used but can add as argument if needed

  my $tr_seq;
  my $floc = get_featureloc(db => $dbh, feature => $tr_id, srctype => $srcfeat_type,);
  #print Dumper($floc);
  my $srcfeat_seq = get_residues($dbh, $floc->{srcfeature_id});

  my $ref = get_exons_of_transcript($dbh, $tr_id, $srcfeat_type, $fl_rank);
  my @allexref = @$ref;
  #print "got exons\n";

  my $counter = 0;
  foreach $exref (@allexref) {
    $counter++;
    undef(my $s);
    undef(my $s1);
    my @ex_ary = @$exref;
    my $ex_fmin = $ex_ary[0];
    my $ex_fmax = $ex_ary[1];
    my $ex_strand = $ex_ary[2];
    my $phase = $ex_ary[3];
    #print "$ex_fmin..$ex_fmax, $ex_strand\n";
    
    if ($counter == 1 and $is_prediction) {
      if ($ex_strand == 1) {
	$s = substr($srcfeat_seq, $ex_fmin + $phase, $ex_fmax - $ex_fmin - $phase);
      } elsif ($ex_strand == -1) {
	$s = substr($srcfeat_seq, $ex_fmin, $ex_fmax - $phase - $ex_fmin);
      }
    } else {
      $s = substr($srcfeat_seq, $ex_fmin, $ex_fmax - $ex_fmin);
    }

    if ($ex_strand == -1) {
      $s1 = reverse_complement($s);
    } else {
      $s1 = $s;
    }
    
    $tr_seq = $tr_seq . $s1;
  }

  $tr_seq = uc($tr_seq);
  return $tr_seq;
}


sub get_cds_seq {
  my $dbh = shift;
  my $cds_id = shift;
  my $cds_srcfeat_type = shift;
  my $cds_seq = '';
  my $fl_rank = shift || 0;
  
  my $floc = get_featureloc(db => $dbh, feature => $cds_id, srctype => $cds_srcfeat_type,);
  my $srcfeat_seq = get_residues($dbh, $floc->{srcfeature_id});
  (my $mrna_id) = get_object_ids(db => $dbh, subject => $cds_id, reltype => 'producedby',);

  my $cds_fmin = $floc->{fmin};
  my $cds_fmax = $floc->{fmax};
  #print "cds fmin: $cds_fmin; cds_fmax: $cds_fmax\n";

  my $ref = get_exons_of_transcript($dbh, $mrna_id, $cds_srcfeat_type, $fl_rank);
  my @allexref = @$ref;

  foreach $exref (@allexref) {
    my @ex_ary = @$exref;
    my $ex_fmin = $ex_ary[0];
    my $ex_fmax = $ex_ary[1];
    my $ex_strand = $ex_ary[2];
    #print "exon fmin: $ex_fmin; fmax: $ex_fmax\n";
    #print "exon len: ", $ex_fmax - $ex_fmin, "\n";
    undef(my $s);
    
    #get the exon coding sequence
    #internal exon
    if ($cds_fmin <= $ex_fmin && $cds_fmax >= $ex_fmax) {
      $s = substr($srcfeat_seq, $ex_fmin, $ex_fmax - $ex_fmin);
    }
    #partial exon
    elsif ($cds_fmin <= $ex_fmin && $cds_fmax <= $ex_fmax && $cds_fmax >= $ex_fmin) {
      $s = substr($srcfeat_seq, $ex_fmin, $cds_fmax - $ex_fmin);
    }
    #partial exon
    elsif ($cds_fmin >= $ex_fmin && $cds_fmin <= $ex_fmax && $cds_fmax >= $ex_fmax) {
      $s = substr($srcfeat_seq, $cds_fmin, $ex_fmax - $cds_fmin);
    }
    #CDS enclosed in this single exon
    elsif ($cds_fmin >= $ex_fmin && $cds_fmax <= $ex_fmax) {
      $s = substr($srcfeat_seq, $cds_fmin, $cds_fmax - $cds_fmin);
    }
    
    if ($ex_strand == -1) {
      $s = reverse $s;
      $s =~ tr/ACGTacgt/TGCAtgca/;
    }
    
    $cds_seq = $cds_seq . $s;
    #print $cds_seq, "\n";
    #print length($cds_seq), "\n";
  }
  
  $cds_seq = uc($cds_seq);
  return $cds_seq;
}



# takes 2 dates both of form yyyy-mm-dd plus optional timestamp of form hour:min:sec
# lack of timestamp assumes 12:00:00.000001 of date provided
# returns true if first date is earlier than 2nd
sub is_earlier_than {
    my $firstdate = shift;
    my $seconddate = shift;

    print "WARNING - you need to supply to dates for an accurate test -- RETURN\n"
      and return unless $firstdate and $seconddate;

    my $date1; my $time1; my $date2; my $time2;
    ($date1, $time1) = split /\s+/, $firstdate if $firstdate;
    ($date2, $time2) = split /\s+/, $seconddate if $seconddate;

    my @date1 = split /-/, $date1;
    my @date2 = split /-/, $date2;

    if ($date1[0] < $date2[0]) {
      return 1;
    } elsif ($date1[1] < $date2[1]) {
      return 1;
    } elsif ($date1[2] < $date2[2]) {
      return 1;
    } else {
      return 1 if $time2 and !$time1;
      return unless $time1 and $time2;
      # if we've got to here need to parse time
      my @time1 = split /:/, $time1;
      my @time2 = split /:/, $time2;
      if ($time1[0] < $time2[0]) {
	return 1;
      } elsif ($time1[1] < $time2[1]) {
	return 1;
      } elsif ($time1[2] < $time2[2]) {
	return 1;
      } else {
	return;
      }
    }
}

sub trim {
    my @s = @_;
    for (@s) {s/^\s+//; s/\s+$//;}
    return wantarray ? @s : $s[0];
}

# returns the number of items in the array if all the items in the array are
# defined or not zero
# NOT CURRENTLY USED IN EvidenceChecker
sub _all_true {
    my @array = @_;
    my $l = @array;
    return $l if (grep$_, @array) == $l;
    0;
}

#returns date in yyyy-mm-dd format
sub get_date {
    my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) 
      = localtime();

    my $year = 1900 + $yearOffset;
    $month++;
    $month = "0$month" if $month < 10;
    $dayOfMonth = "0$dayOfMonth" if $dayOfMonth < 10;

    return "$year-$month-$dayOfMonth";
}

#creates a string of feature_ids to use in select statements to limit
#queries to only specified features - only works in get_gene_info query?
sub create_id_string {
  my @id_list = @{+shift};
  $id_string = "and f3.feature_id IN (";
  for (0..$#id_list) {
    $id_string .= "$id_list[$_] ";
    if ($_ == $#id_list) {
      $id_string .= ")";
    } else {
      $id_string .= ", ";
    }
  }
  $id_string;
}
    
# creates a hashref with keys = dereferenced array_ref and 
# values = supplied scalars
# arg1 - arrayref of keys
# arg2..n = scalar values
# return reference to created hash
sub _create_hashref {
    my %hash;
    my @keys = @{+shift};
    
    foreach $key (@keys) {
	$hash{$key} = shift;
    }
    \%hash;
}

# checks if 2 arrays are identical
#
# args - array refs to 2 arrays to compare
# return true if identical
sub arrays_r_identical {
     my @a1 = @{+shift};
     my @a2 = @{+shift};

     return unless scalar @a1 == scalar @a2;

     for(0..$#a1) {
	  return unless $a1[$_] eq $a2[$_];
     }

     return 1;
}


# sorts genes into position based on arm, strand and fmin
# we do care about strand
sub by_arm_strand_min {
    $a->{src} <=> $b->{src}
      or
    $a->{strand} <=> $b->{strand}
	  or
    $a->{min} <=> $b->{min};
}

# sorts genes into position based on arm and fmin
# we don't care about strand
sub by_arm_min {
    $a->{src} <=> $b->{src}
      or
    $a->{min} <=> $b->{min};
}

sub by_srcfeature_id_fmin {
  $a->{srcfeature_id} <=> $b->{srcfeature_id}
      or
    $a->{fmin} <=> $b->{fmin};
}


#sorting functions
#by fmin or fmax
sub by_min {$a->{min}<=>$b->{min}}
sub by_fmin {$a->{fmin}<=>$b->{fmin}}
sub by_max {$b->{max}<=>$a->{max}}
#from 5' to 3' regardless of strand
sub by_5to3 {($a->{min} * $a->{strand}) <=> ($b->{min} * $b->{strand})} 


sub reverse_complement {
  my @revcomp = map { (my $rev = reverse $_) =~ tr /ACGTacgt/TGCAtgca/; $rev } @_;
  return $revcomp[0] if scalar @revcomp == 1;
  return @revcomp;
}

sub is_protein {
    my $seq = shift;

    return if $seq =~ /[^acdefghiklmnpqrstvwy]/i;
    return 1;
}

sub is_protein_ambig {
    my $seq = shift;

    return if $seq =~ /[^abcdefghiklmnpqrstuvwxyz\*]/i;
    return 1;
}


sub is_rna
{
  my $seq = shift;
  return if $seq =~ /[^acgu]/i;
  return 1;
}

sub is_dna
{
    my $seq = shift;

    return if $seq =~ /[^acgt]/i;
    return 1;
}


sub is_dna_ambig
{
  my $seq = shift;

  return if $seq =~ /[^acgtn]/i;
  return 1;
}



sub is_nucleic
{
    my $seq = shift;

    return if $seq =~ /[^acgtnmrwsykvhdbx]/i;
    return 1;
}

sub from_gff_string {
    my ( $feat, $string) = @_;
    #print $string;
    chomp($string);

    # according to the nascent GFF3 spec, it should be space
    # separated elements, spaces anywhere else should be escaped.

    my ($seqname, $source, $primary, $start, $end,
   $score, $strand, $frame, $groups) = split(/\t/, $string);

    if ( ! defined $frame ) {
      print "[$string] does not look like GFF3 to me\n";
    }
    $feat->seq_id($seqname);
    $feat->source_tag($source);
    $feat->primary_tag($primary);
    $feat->start($start);
    $feat->end($end);
    $feat->frame($frame);
    if ( $score eq '.' ) {
      #$feat->score(undef);
    } else {
      $feat->score($score);
    }
    if ( $strand eq '-' ) { $feat->strand(-1); }
    if ( $strand eq '+' ) { $feat->strand(1); }
    if ( $strand eq '.' ) { $feat->strand(0); }
    my @groups = split(/\s*;\s*/, $groups);

    for my $group (@groups) {
   my ($tag,$value) = split /=/,$group;
   $tag             = unescape($tag);
   my @values       = map {unescape($_)} split /,/,$value;
   for my $v ( @values ) {  $feat->add_tag_value($tag,$v); }
    }
}
sub unescape {
  my $v = shift;
  #$v =~ tr/+/ /;
  $v =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $v;
}

sub get_highest_caller_line_number {
    my $level = 0;
    my $line_number;
    while (my @caller_info = caller($level++)) {
        $line_number = $caller_info[2];
    }
    return $line_number;
}

sub get_formatted_time {
  my $current_time = localtime();
  my $formatted_time = $current_time->strftime("%Y-%m-%d %H:%M:%S");
  return $formatted_time;
}

sub print_log {
  my $text = shift;
  my $time = get_formatted_time();
  my $line_number = get_highest_caller_line_number();
  print "$time : Line No $line_number : $text\n";
}

1;
