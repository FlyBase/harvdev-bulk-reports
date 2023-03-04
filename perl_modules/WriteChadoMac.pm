#!/usr/bin/perl -w
# package FlyBase::WriteChado;

# use 5.008005;
use strict;
use XML::DOM;

=head1 NAME

  WriteChadoMac.pm - A module to write chado xml elements

  Updated Version - can be used to produce macroized chado-xml

=head1 SYNOPSIS

 use XML::DOM;
 use WriteChadoMac;
 use PrettyPrintDom;
 $doc = new XML::DOM::Document;
 $feat_el = create_ch_feature(doc => $doc,
                              uniquename => 'Hoppy',
                              genus => 'Bufo'
                              species => 'marinus'
                              type => 'invader'
                              with_id => 1,
                              );

 pretty_print($feat_el,\*STDOUT);

 This module can be used to produce either verbose or macroized chado-xml

 NOTE: that this is not backward compatible with previous version that only 
       produced verbose chado-xml 

 Many of the elements are keyed by column name.
 HOWEVER, in some cases where a column is an _id column referencing fields from a different table
 then if the parameter is specified as column name (e.g. type_id) then an id referencing a previously
 defined element or an XML::DOM element itself is expected as the argument.  In a subset of these cases 
 there is a parameter of the same name as the _id column lacking the _id suffix (e.g. type).  These expect
 a string which will be converted into an element as appropriate and specified which is then added
 to the parent element.

 Check the method descriptions for allowed parameters.

 When producing macroized xml the caller is responsible for ensuring that a macro id 
 is assigned for an element that will be used later and must keep track of these ids.
 
 WARNING - for the most part garbage in means garbage out as there is not tons of
           error checking implemented

 FUNCTIONS EXIST TO CREATE CHADO XML FOR THE FOLLOWING 

 analysis and analysis_id
 analysisprop
 analysisfeature
 analysisgrp
 analysisgrpmember
 cell_line
 cell_lineprop
 cell_lineprop_pub
 cell_line_cvterm
 cell_line_cvtermprop
 cell_line_dbxref
 cell_line_feature
 cell_line_library
 cell_line_libraryprop
 cell_line_pub
 cell_line_relationship
 cell_line_synonym
 contact and contact_id
 cv and cv_id
 cvterm
 cvtermprop
 cvterm_relationship
 db and db_id
 dbxref and dbxref_id
 dbxrefprop
 environment and environment_id
 environment_cvterm
 expression and expression_id
 expression_cvterm
 expression_cvtermprop
 expression_pub
 expressionprop
 feature
 feature_cvterm
 feature_cvtermprop
 feature_dbxref
 feature_expression
 feature_expressionprop
 feature_interaction
 feature_interactionprop
 feature_interaction_pub
 feature_genotype
 feature_grpmember
 feature_grpmember_pub
 feature_pub
 feature_pubprop
 feature_relationship (can be either subject or object)
 feature_relationshipprop
 feature_relationshipprop_pub
 feature_relationship_pub
 feature_synonym
 featureloc
 featureloc_pub
 featureprop
 featureprop_pub
 genotype and genotype_id
 grp
 grp_cvterm
 grp_dbxref
 grpprop
 grpprop_pub
 grp_pub
 grp_pubprop
 grp_relationship
 grp_relationship_pub
 grp_relationshipprop
 grp_synonym
 grpmember
 grpmember_pub
 grpmemberprop
 grpmemberprop_pub
 grpmember_cvterm
 humanhealth
 humanhealth_cvterm
 humanhealth_dbxref
 humanhealth_dbxrefprop
 humanhealth_dbxrefprop_pub
 humanhealth_feature
 humanhealth_featureprop
 humanhealth_relationship
 humanhealth_relationship_pub
 humanhealth_phenotype
 humanhealth_phenotypeprop 
 humanhealth_pub
 humanhealth_pubprop
 humanhealth_synonym
 humanhealthprop
 humanhealthprop_pub
 interaction and interaction_id
 interactionprop
 interactionprop_pub
 interaction_cell_line
 interaction_cvterm
 interaction_cvtermprop
 interaction_expression
 interaction_expressionprop
 interaction_pub
 library and library_id
 library_cvterm
 library_cvtermprop
 library_dbxref
 library_dbxrefprop
 libraryprop
 libraryprop_pub
 library_expression
 library_expressionprop
 library_feature
 library_featureprop
 library_grpmember
 library_humanhealth
 library_humanhealthprop
 library_interaction
 library_pub
 library_relationship
 library_relationship_pub
 library_strain
 library_strainprop
 library_synonym
 organism and organism_id
 organism_cvterm
 organism_cvtermprop
 organism_dbxref
 organism_grpmember
 organism_library
 organism_pub
 organismprop
 organismprop_pub
 phendesc
 phenotype and phenotype_id
 phenotype_comparison
 phenotype_comparison_cvterm
 phenotype_cvterm
 phenstatement
 pub and pub_id
 pubauthor
 pubprop
 pub_dbxref
 pub_relationship
 strain
 strain_cvterm
 strain_dbxref
 strain_feature
 strain_featureprop
 strain_relationship
 strain_relationship_pub
 strain_phenotype
 strain_phenotypeprop 
 strain_pub
 strain_synonym
 strainprop
 strainprop_pub
 synonym and synonym_id
 generic prop


=head1 DESCRIPTION

=head2 Methods

=over 12

=item C<create_ch_analysis>

 CREATE analysis or analysis_id element
 params 
 doc - XML::DOM::Document optional - required
 program - required string
 programversion - required string (usually a number) NOTE: can add default 1.0?
 sourcename - optional string NOTE: this is part of the unique key and while it can be it usually shouldn't be null
 name - optional string
 description - optional string
 algorithm - optional string
 sourceversion - optional string
 sourceuri -optional string
 timeexecuted - optional string value like '1999-01-08 04:05:06' default will be whenever data is added
                NOTE: not sure on xort-postgres interaction regarding invalid timestamp formats
 macro_id - string optional if provide then add an ID attribute to the top level element of provided value
 with_id - optional if true will create analysis_id at top level

=item C<create_ch_analysisfeature>

 CREATE analysisfeature element
 params 
 doc - XML::DOM::Document optional - required

 Parameters to make a feature element must either pass a feature_id or the other necessary bits
 feature_id - macro feature id or XML::DOM feature element
 uniquename - string
 organism_id - macro organism id or XML::DOM organism element
 genus - string
 species - string
 type_id -  macro id for feature type or XML::DOM cvterm element for type
 type - string valid SO feature type

 Parameters to make an analysis element must either pass analysis_id or required bits
 analysis_id - macro analysis id or XML::DOM analysis element
 program - string
 programversion - string
 sourcename - string

 Here are the optional bits that can be added to the analysisfeature
 rawscore - number (double) 
 normscore - number (double)  
 significance - number (double)
 identity - number (double)    


=item C<create_ch_analysisprop>

 CREATE analysisprop element
 params
 doc - XML::DOM::Document required
 analysis_id - optional macro id for an analysis or XML::DOM analysis element for standalone analysisprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for a analysis property type or XML::DOM cvterm element required
 type - string from analysis property type cv
        Note: will default to making a featureprop from 'analysis property type' cv unless cvname is provided
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_analysisgrp>

 CREATE analysisgrp element
 params 
 doc - XML::DOM::Document optional - required

 Parameters to make a grp element must either pass a grp_id or the other necessary bits
 grp_id - macro grp id or XML::DOM grp element
 uniquename - string
 type_id -  macro id for grp type or XML::DOM cvterm element for type
 type - string valid SO grp type

 Parameters to make an analysis element must either pass analysis_id or required bits
 analysis_id - macro analysis id or XML::DOM analysis element
 program - string
 programversion - string
 sourcename - string

 Here are the optional bits that can be added to the analysisgrp
 rawscore - number (double) 
 normscore - number (double)  
 significance - number (double)
 identity - number (double)    

=item C<create_ch_analysisgrpmember>

 CREATE analysisgrpmember element
 params
 doc - XML::DOM::Document required
 analysis_id  optional macro id for a analysis or XML::DOM analysis element for standalone analysis_grpmember
 grpmember_id - optional macro id for a grpmember or XML::DOM grpmember element for standalone analysisgrpmember

=item C<create_ch_cell_line>

 CREATE cell_line element
 params
 doc - XML::DOM::Document required
 uniquename - string required
 organism_id - organism macro id or XML::DOM organism element required if no genus and species
 genus - string required if no organism
 species - string required if no organism
 name - string optional
 with_id - boolean optional if 1 then cell_line_id element is returned
 no_lookup - boolean option if 1 then default op="lookup" attribute will not be added to element

=item C<create_ch_cell_lineprop>

 CREATE cell_lineprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - string from  cell_lineprop type cv or XML::DOM cvterm element required 
        Note: will default to making a featureprop from 'cell_lineprop type' cv unless 
              cvname is provided
 cvname - string optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_cell_lineprop_pub>

 CREATE cell_lineprop_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of cell_lineprop_pub
      or just appending the pub element to cell_lineprop_pub if that is passed
 params
 doc - XML::DOM::Document required
 pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_cell_line_cvterm>

 CREATE cell_line_cvterm element
 params
 doc - XML::DOM::Document required
 cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
 name - string required unless cvterm_id provided Note: a cvterm has a lookup by default cannot make a new cvterm 
                                                        with this method
 cv_id - macro id for cv or XML::DOM cv element required if name and not cv
 cv - string for name of cv required if name and not cv_id
 pub_id - macro id for pub or XML::DOM pub element required unless pub
 pub - string = pub uniquename Note: as pub has lookup option by default can't make a new pub using this param
 rank - int optional default = 0

=item C<create_ch_cell_line_cvtermprop>

 CREATE cell_line_cvtermprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - string from  cell_line_cvtermprop type cv or XML::DOM cvterm element required 
        Note: will default to making a cell_lineprop from 'cell_line_cvtermprop type' cv unless 
              cvname is provided
 cvname - string optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_cell_line_dbxref>

 CREATE cell_line_dbxref element
 params
 doc - XML::DOM::Document required
 cell_line_id - macro cell_line id or XML::DOM cell_line element optionaal to create freestanding cell_line_dbxref
 dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
 accession - string required unless dbxref_id provided
 db_id - macro db id or XML::DOM db element required unless dbxref_id provided
 db - string name of db
 version - string optional
 description - string optional
 is_current - string 't' or 'f' boolean default = 't' so don't pass unless
              this shoud be changed

=item C<create_ch_cell_line_feature>

 CREATE cell_line_feature element
 params
 doc - XML::DOM::Document required
 organism_id - macro id for organism or XML::DOM organism element
 genus - string
 species - string
 NOTE: you can use the generic paramaters in the following cases:
       1.  you are only building either a cell_line or feature element and not both
       2.  or both cell_line and feature have the same organism
       otherwise use the prefixed parameters
 WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
 cell_line_id - macro id for cell_line or XML::DOM cell_line element
 cell_uniquename - string cell_line uniquename
 cell_organism_id - macro id for organism or XML::DOM organism element to link to cell_line
 cell_genus
 cell_species
 feature_id - macro id for feature or XML::DOM feature element
 feat_uniquename
 feat_organism_id - macro id for organism or XML::DOM organism element to link to feature
 feat_genus
 feat_species
 feat_type_id
 feat_type
 pub_id - macro pub id or XML::DOM pub element - required unless puname provided
 pub - string uniquename for pub (note will have lookup so can't create new pub here)

=item C<create_ch_cell_line_library>

 CREATE cell_line_library element
 params
 doc - XML::DOM::Document required
 organism_id - macro id for organism or XML::DOM organism element
 genus - string
 species - string
 NOTE: you can use the generic paramaters in the following cases:
       1.  you are only building either a cell_line or library element and not both
       2.  or both cell_line and library have the same organism
       otherwise use the prefixed parameters
 WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
 cell_line_id - macro id for cell_line or XML::DOM cell_line element
 cell_uniquename - string cell_line uniquename
 cell_organism_id - macro id for organism or XML::DOM organism element to link to cell_line
 cell_genus
 cell_species
 library_id - macro id for library or XML::DOM library element
 lib_uniquename
 lib_organism_id - macro id for organism or XML::DOM organism element to link to library
 lib_genus
 lib_species
 lib_type_id
 lib_type
 pub_id - macro pub id or XML::DOM pub element - required unless puname provided
 pub - string uniquename for pub (note will have lookup so can't create new pub here)

=item C <create_ch_cell_line_libraryprop>

 CREATE cell_line_libraryprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for cell_line_libraryprop type or XML::DOM cvterm element required
 type -  string from  cell_line_libraryprop type cv 
        Note: will default to making a cell_line_libraryprop from 'cell_line_libraryprop type' cv unless cvname is provided
 cvname - string (probably want to pass 'cell_line_libraryprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0


=item C<create_ch_cell_line_pub>

 CREATE cell_line_pub element
 Note that this is just calling create_ch_pub  
      and adding returned pub_id element as a child of cell_line_pub
      or just appending the pub element to cell_line_pub if that is passed
 params
 doc - XML::DOM::Document required
 pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
        creating a new pub (i.e. not null value but not part of unique key
 type - string for type from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_cell_line_relationship>

 CREATE cell_line_relationship element
 NOTE: this will create either a subject_id or object_id cell_line_relationship element
 but you have to attach this to the related cell_line elsewhere
 params
 doc - XML::DOM::Document required
 object _id - macro id for object cell_line or XML::DOM cell_line element
 subject_id - macro id for subject cell_line or XML::DOM cell_line element
 NOTE you can pass one or both of the above parameters with the following rules:
 if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other cell_line
 if both are passed then is_object, is_subject and any parameters to create a cell_line are ignored
 is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
 is_subject - boolean 't'          this flag indicates if the cell_line info provided should be 
                                   added in as subject or object cell_line
 type_id - macro id for relationship type or XML::DOM cvterm element (Note: currently all is_relationship = '0'
 type - string for relationship type note: with this param  type will be assigned to relationship_type cv
 rank - integer optional with default = 0
 cell_line_id - macro_id for a cell_line or XML::DOM cell_line element required unless minimal cell_line bits provided
 uniquename - string required unless cell_line provided
 organism_id - macro id for organism or XML::DOM organism element required unless cell_line or (genus & species) provided
 genus - string required unless cell_line or organism provided
 species - string required unless cell_line or organism provided

 Alias: create_ch_clr
 Alias: create_ch_c_l_r

=item C<create_ch_cell_line_synonym>

 CREATE cell_line_synonym element
 params
 doc - XML::DOM::Document required
 synonym_id - XML::DOM synonym element required unless name and type provided
 name - string required unless synonym_id element provided
 type_id - macro id for synonym type or XML::DOM cvterm element
                 required unless a synonym element provided
 type - string = name from the 'synonym type' cv
 pub_id macro id for a pub or a XML::DOM pub element required
 pub - a pub uniquename (i.e. FBrf)
 synonym_sgml - string optional but if not provided then synonym_sgml = name
               - do not provide if a synonym element is provided
 is_current - optional string = 'f' or 't' default is 't' so don't provide this param
               unless you know you want to change the value
 is_internal - optional string = 't' default is 'f' so don't provide this param
               unless you know you want to change the value

=item C<create_ch_contact>

 CREATE contact or contact_id element
 params
 doc - XML::DOM::Document optional - required
 name - string required
 description - string optional
 macro_id - optional string to specify as ID attribute for contact
 with_id - boolean optional if 1 then contact_id element is returned

=item C<create_ch_cv>

 CREATE cv or cv_id element
 params
 doc - XML::DOM::Document required
 name - string required
 definition - string optional
 macro_id - optional string to specify as ID attribute for cv
 with_id - boolean optional if 1 then cv_id element is returned

=item C<create_ch_cvterm>

 CREATE cvterm element
 params
 doc - XML::DOM::Document required
 name - string
 cv_id - macro id for cv or XML::DOM cv element
 cv - string = cvname
 definition - string optional
 dbxref_id - macro id for dbxref XML::DOM dbxref element
 is_obsolete - boolean optional default = 0
 is_relationshiptype - boolean optional default = 0
 macro_id - optional string to specify as ID attribute for cvterm
 no_lookup - boolean optional
 note that we don't have a with_id parameter because either it will be freestanding 
 term or will have another type of id (e.g. type_id)

 note: there are 2 unique keys on the cvterm table (name, cv_id, is_obsolete) and (dbxref_id)
 this method requires that all the info for at least one of the unique keys is present
 it is up to the to make sure that the right key is used upon loading

=item C<create_ch_cvtermprop>

 CREATE cvtermprop element
 params
 doc - XML::DOM::Document required
 cvterm_id - optional macro id for a cvterm or XML::DOM cvterm element for standalone cvtermprop
 value - string - not strictly required and in some cases this value is null in chado
 type_id - macro id for cvterm property type or XML::DOM cvterm element required
 type - string from cvterm_property_type cv
        Note: will default to making a cvtermprop from above cv unless cvname is provided
 cvname - string optional 
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_cvterm_relationship>

 CREATE cvterm_relationship element
 NOTE: this can now create a free standing cvterm relationship if you pass subject_id or object_id
 params
 doc - XML::DOM::Document required
 object _id - macro id for object feature or XML::DOM feature element
 subject_id - macro id for subject feature or XML::DOM feature element
 is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
 is_subject - boolean 't'          this flag indicates if the cvterm info provided should be 
                                   added in as subject or object cvterm
 rtype_id -  macro for relationship type or XML::DOM cvterm element 
 rtype - string for relationship type note: if relationship name is given will be assigned to relationship_type cv
 (Note: currently all is_relationship = '0')
 cvterm_id - macro id for cvterm or XML::DOM cvterm element required unless name and cv info provided
 name - string
 cv_id - macro id for a cv or XML::DOM cv element
 cv - cv name string required if name and not cv_id
 dbxref_id - macro id for dbxref or XML::DOM dbxref element
 macro_id - optional string to add as ID attribute value to cvterm

 Alias: create_ch_cr

=item C<create_ch_db>

 CREATE db or db_id element
 params
 doc - XML::DOM::Document required
 name - string required
 contact_id - macro id for contact or XML::DOM contact element optional
 contact - string = contact name
 description - string optional
 urlprefix - string optional
 url - string optional
 macro_id - optional string to add as ID attribute value to db
 with_id - boolean optional if 1 then db_id element is returned

=item C<create_ch_dbxref>

 CREATE dbxref or dbxref_id element
 params
 doc - XML::DOM::Document required
 accession - string required
 db_id - macro id for db or XML::DOM db element required unless db
 db - string = db name required unless db_id
 version - string optional will default to ''
 description - string optional
 macro_id - optional string to add as ID attribute value to dbxref
 with_id - boolean optional if 1 then dbxref_id element is returned
 no_lookup - boolean optional

=item C<create_ch_dbxrefprop>

 CREATE dbxrefprop element
 params
 doc - XML::DOM::Document required
 dbxref_id - optional macro id for a dbxref or XML::DOM dbxref element for standalone dbxrefprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for property type or XML::DOM cvterm element required
 type - string from dbxrefprop type cv 
        Note: will default to making a featureprop from 'dbxrefprop type' cv unless cvname is provided
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0


=item C<create_ch_environment>

 CREATE environment or environment_id element
 params
 doc - XML::DOM::Document optional - required
 uniquename - string required
 description - string optional
 with_id - boolean optional if 1 then expression_id element is returned
 macro_id - string to specify as id attribute of this element for later use

=item C<create_ch_environment_cvterm>

 CREATE environment_cvterm element
 params
 doc - XML::DOM::Document optional - required
 environment_id - macro id for environment of XML::DOM environment element
 uniquename - environment uniquename
 NOTE: you need to pass environment bits if attaching to existing cvterm element or 
       creating a freestanding environment_cvterm
 cvterm_id -  macro id for cvterm of XML::DOM cvterm element
 name - cvterm name
 cv_id - macro id for a CV or XML::DOM cv element
 cv - name of a cv
 is_obsolete - optional param for cvterm
 NOTE: you need to pass cvterm bits if attaching to existing environment element or 
       creating a freestanding environment_cvterm

=item C<create_ch_expression>

 CREATE expression or expression_id element
 params
 doc - XML::DOM::Document
 uniquename - string required
 description - string optional
 md5checksum - char(32) optional
 macro_id - optional string to add as ID attribute value to expression
 with_id - boolean optional if 1 then expression_id element is returned

=item C<create_ch_expression_cvterm>

 CREATE expression_cvterm 
 params
 doc - XML::DOM::Document - required
 expression_id - OPTIONAL macro expression id or XML::DOM expression element to create freestanding expression_cvterm
 cvterm_id - macro id for a cvterm or XML::DOM cvterm element - required unless name and cv params
 name - string name for cvterm required unless cvterm_id
 cv_id - macro id for cv or XML::DOM cv element required unless cvterm_id or cv provided
 cv - string = cvname required unless cvterm_id or cv_id provided
 cvterm_type_id - macro id for expression slot cvterm or XML::DOM cvterm element - required unless type
 cvterm_type - string from the expression slots cv
 rank - integer with default of 0 - this plus type_id used for ordering items 

=item C<create_ch_expression_cvtermprop>

 CREATE expression_cvtermprop
 params
 doc - XML::DOM::Document - require
 expression_cvterm_id - optional macro id for a expression_cvterm or XML::DOM expression_cvterm element 
                        for standalone expression_cvtermprop
 value - string
 type_id - macro id for a expression_cvtermprop type cvterm or XML::DOM cvterm element required
 type - string from expression_cvtermprop type cv
     Note: will default to making a cvtermprop from 'expression_cvtermprop type' cv unless cvname is provided
 cvname - string optional but see above for the default 
 rank - integer with a default of zero so don't use unless you want a rank other than 0


=item C<create_ch_expression_pub>

 CREATE expression_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of expression_pub
      or just appending the pub element to expression_pub if that is passed 
 params
 doc - XML::DOM::Document required
 expression_id - OPTIONAL macro expression id or XML::DOM expression element to create freestanding expression_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc (optional expression_id)  as only params
 uniquename - string for pub uniquename required unless pub_id
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
           creating a new pub (i.e. not null value but not part of unique key
 type - string from pub type cv same requirement as type_id
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_expressionprop>

 CREATE expressionprop element
 params
 doc - XML::DOM::Document required
 expression_id - optional macro id for a expression or XML::DOM expression element for standalone expressionprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for property type or XML::DOM cvterm element required
 type - string from expressionprop type cv 
        Note: will default to making a featureprop from 'expressionprop type' cv unless cvname is provided
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0


=item C<create_ch_feature>

 CREATE feature element
 params
 doc - XML::DOM::Document required
 uniquename - string required
 type_id - macro id for cvterm or XML::DOM cvterm element required 
 type - string for a cvterm name Note: will default to using SO cv unless a cvname is provided
 cvname - string optional to specify a cv other than SO for the type_id
          do not use if providing a cvterm element to type_id
 organism_id - organism macro id or XML::DOM organism element required if no genus and species
 genus - string required if no organism
 species - string required if no organism
 dbxref_id - dbxref macro id or XML::DOM dbxref element optional
 name - string optional
 residue - string optional
 seqlen - integer optional (if seqlen = 0 pass as string)
 md5checksum - string optional
 is_analysis - boolean 't' or 'f' default = 'f' optional
 is_obsolete - boolean 't' or 'f' default = 'f' optional
 macro_id - optional string to add as ID attribute value to feature
 with_id - boolean optional if 1 then feature_id element is returned
 no_lookup - boolean optional

=item C<create_ch_feature_cvterm>

 CREATE feature_cvterm element
 params
 doc - XML::DOM::Document required
 feature_id - OPTIONAL macro feature id or XML::DOM feature element to create freestanding feature_cvterm
 cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
 name - string required unless cvterm_id provided
 cv_id - macro id for cv or XML::DOM cv element required unless cvterm_id provided
 cv - string = cvname
 pub_id - macro id for a pub or XML::DOM pub element required unless pub 
 pub - string pub uniquename Note: cannot create a new pub if uniquename provided
 is_not - optional boolean 't' or 'f' with default = 'f' so don't pass unless you know you want to change

=item C<create_ch_feature_cvtermprop>

 CREATE feature_cvtermprop element
 params
 doc - XML::DOM::Document required
 feature_cvterm_id - optional macro id for a feature_cvterm or XML::DOM feature_cvterm element 
                     for standalone feature_cvtermprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for feature_cvtermprop type or XML::DOM cvterm element required
 type -  string from  feature_cvtermprop type cv 
        Note: will default to making a featureprop from 'property type' cv unless cvname is provided
 cvname - string (probably want to pass 'feature_cvtermprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_feature_dbxref>

 CREATE feature_dbxref element
 params
 doc - XML::DOM::Document required
 feature_id - OPTIONAL macro feature id or XML::DOM feature element to create freestanding feature_dbxref
 dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
 accession - string required unless dbxref_id provided
 db_id - macro db id or XML::DOM db element required if accession and not db provided
 db - string = db name required if accession and not db_id provided
 version - string optional
 description - string optional
 is_current - string 't' or 'f' boolean default = 't' so don't pass unless
              this should be changed

=item C<create_ch_feature_expression>

 CREATE feature_expression element
 params
 doc - XML::DOM::Document required
 feature_id - OPTIONAL macro feature id or XML::DOM feature element to create freestanding feature_expression 
 expression_id - macro expression id or XML::DOM expression element - required unless uniquename provided
 uniquename - string required unless expression_id
 pub_id - macro pub id or XML::DOM pub element - required unless puname provided
 pub - string uniquename for pub (note will have lookup so can't create new pub here)

=item C<create_ch_feature_expressionprop>

 CREATE feature_expressionprop element
 params
 doc - XML::DOM::Document required
 feature_expression_id - optional macro id for a feature_expression or XML::DOM feature_expression element 
                         for standalone feature_expressionprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for feature_expressionprop type or XML::DOM cvterm element required
 type -  string from  expression_cvtermprop type cv 
        Note: will default to making a featureprop from 'expression_cvtermprop type' cv unless cvname is provided
 cvname - string (probably want to pass 'feature_expressionprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_feature_genotype>

 CREATE feature_genotype or feature_genotype_id element
 params
 doc - XML::DOM::Document - required
 feature_id - OPTIONAL macro feature id or XML::DOM feature element to create freestanding feature_genotype
 genotype_id - macro id for genotype or XML::DOM genotype element required unless uniquename
 uniquename - string = genotype uniquename
 chromosome_id - macro id for chromosome feature or XML::DOM feature element required
 cvterm_id - macro id for a cvterm or XML::DOM cvterm element required 
 rank - integer optional with default = 0
 cgroup - integer optional with default = 0

=item C<create_ch_feature_grpmember>

 CREATE feature_grpmember element
 params
 doc - XML::DOM::Document required
 feature_id  optional macro id for a feature or XML::DOM feature element for standalone feature_grpmember
 grpmember_id - optional macro id for a grpmember or XML::DOM grpmember element for standalone feature_grpmember

=item C<create_ch_feature_grpmember_pub>

 CREATE feature_grpmember_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of feature_grpmember_pub
      or just appending the pub element to feature_grpmember_pub if that is passed
 params
 doc - XML::DOM::Document required
 feature_grpmember_id - OPTIONAL macro feature_grpmember id or XML::DOM feature_grpmember element to create freestanding feature_grpmember_pub
 pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_feature_humanhealth_dbxref>

 CREATE feature_humanhealth_dbxref element
 params
 doc - XML::DOM::Document required
 feature_id - OPTIONAL macro feature id or XML::DOM feature element to create freestanding feature_humanhealth_dbxref 
 humanhealth_dbxref_id - macro humanhealth_dbxref_id or XML::DOM humanhealth_dbxref element - required
 pub_id - macro pub id or XML::DOM pub element - required unless puname provided
 pub - string uniquename for pub (note will have lookup so can't create new pub here)

=item C<create_ch_feature_interaction>

 CREATE a feature_interaction element
 params
 doc - XML::DOM::Document required
 feature_id - OPTIONAL macro feature id or XML::DOM feature element to create freestanding feature_expression 
 interaction_id - macro interaction_id or XML::DOM interaction element - required unless uniquename and type info provided
 uniquename - string required unless interaction_id
 type_id - macro type_id or cvterm  element required unless interaction_id or type
 type - string required unless interaction_id or type_id
 cvname - string optional
 role_id - macro role_id or cvterm element unless role
 role - string term from 'PSI-MI' cv
 rank - int optional default = 0

=item C<create_ch_feature_interactionprop>

 CREATE a feature_interactionprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - string from interaction property type cv or XML::DOM cvterm element required
        Note: will default to making a from 'interaction property type' cv unless cvname is provided
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_feature_interaction_pub>

 CREATE feature_interaction_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of feature_interaction_pub
      or just appending the pub element to feature_interaction_pub if that is passed
 params
 doc - XML::DOM::Document required
 feature_interaction_id - optional feature_interaction element or macro_id 
          to create freestanding ele
 pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_feature_pub>

 CREATE feature_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of feature_pub
      or just appending the pub element to feature_pub if that is passed 
 params
 doc - XML::DOM::Document required
 feature_id - OPTIONAL macro feature id or XML::DOM feature element to create freestanding feature_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc (optional feature_id)  as only params
 uniquename - string for pub uniquename required unless pub_id
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
           creating a new pub (i.e. not null value but not part of unique key
 type - string from pub type cv same requirement as type_id
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_feature_pubprop>

 CREATE feature_pubprop element
 params
 doc - XML::DOM::Document required
 feature_pub_id - optional macro id for a feature_pub or XML::DOM feature_pub element for standalone feature_pubprop
 value - string - not strictly required
 type_id - macro id for feature_pubprop type or XML::DOM cvterm element required 
 type - string from feature_pubprop type cv
        Note: will default to making a pubprop from above cv unless cvname is provided
 cvname - string optional 
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_fr>

 CREATE feature_relationship element
 NOTE: this can now create a free standing feature relationship if you pass subject_id or object_id
 params
 doc - XML::DOM::Document required
 object _id - macro id for object feature or XML::DOM feature element
 subject_id - macro id for subject feature or XML::DOM feature element
 NOTE you can pass one or both of the above parameters with the following rules:
    if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other feature
    if both are passed then is_object, is_subject and any parameters to create a feature are ignored
 is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
 is_subject - boolean 't'          this flag indicates if the feature info provided should be 
                                   added in as subject or object feature
 rtype_id - macro id for cvterm or XML::DOM cvterm element
 rtype - string for relationship type note: if relationship name is given will be assigned to 'relationship type' cv
 rank - integer optional with default = 0
 feature_id - macro id for feature or XML::DOM feature element required unless minimal feature bits provided
 uniquename - string required unless feature_id provided
 organism_id - macro id for an organism or XML::DOM organism element required unless feature_id or (genus & species) provided
 genus - string required unless feature_id or organism_id provided
 species - string required unless feature_id or organism_id provided
 ftype_id -  macro id for cvterm or XML::DOM cvterm element required unless feature provided
 ftype - string = name of feature type 

 Alias: create_ch_feature_relationship
 Alias: create_ch_f_r

=item C<create_ch_frprop>

 CREATE feature_relationshipprop element
 params
 doc - XML::DOM::Document required
 feature_relationship_id - optional macro id for a feature_relationship or XML::DOM feature_relationship element 
                           for standalone feature_relationshipprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for a feature_relationshipprop type or XML::DOM cvterm element required
 type - string from feature_relationshipprop type cv 
        Note: will default to making a featureprop from above cv unless cvname is provided
 cvname - string optional 
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

 Alias: create_ch_fr_prop is a synonym for backward compatibility
 Alias: create_ch_feature_relationshipprop
 Alias: create_ch_f_rprop

=item C<create_ch_frprop_pub>

 CREATE feature_relationshipprop_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of feature_relationshipprop_pub
      or just appending the pub element to feature_relationshipprop_pub if that is passed
 params
 doc - XML::DOM::Document required
 feature_relationshipprop_id -optional feature_relationshipprop XML::DOM element or macro id
 pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

 Alias: create_ch_feature_relationshipprop_pub
 Alias: create_ch_f_rprop_pub

=item C<create_ch_fr_pub>

 CREATE feature_relationship_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of feature_relationship_pub
      or just appending the pub element to feature_relationship_pub if that is passed
 params
 doc - XML::DOM::Document required
 feature_relationship_id -optional feature_relationship XML::DOM element or macro id
 pub_id - XML::DOM pub element - if this is used then pass this and doc (and optional fr_id) as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

 Alias: create_ch_feature_relationship_pub

=item C<create_ch_feature_synonym>

 CREATE feature_synonym element
 params
 doc - XML::DOM::Document required
 feature_id - OPTIONAL macro feature id or XML::DOM feature element to create freestanding feature_synonym
 synonym_id - macro id for synonym or XML::DOM synonym element required unless name and type provided
 name - string required unless synonym_id element provided
 type_id - macro id for cvterm or XML::DOM cvterm element required if name and not type
 type - string = name from the 'synonym type' cv required if name and type_id
 pub_id - macro id for a pub or a XML::DOM pub element required
 pub - a pub uniquename (i.e. FBrf)
 synonym_sgml - string optional but if not provided then synonym_sgml = name
               - do not provide if a synonym_id element is provided
 is_current - optional string = 'f' or 't' default is 't' so don't provide this param
               unless you know you want to change the value
 is_internal - optional string = 't' default is 'f' so don't provide this param
               unless you know you want to change the value

=item C<create_ch_featureloc>

 CREATE featureloc element
 params
 none of these parameters are strictly required and some warning is done
 but if you misbehave you could set up some funky featurelocs?
 srcfeature_id macro id for a feature or a XML::DOM feature element
 fmin - integer (NOTE: if fmin = 0 you can pass as string to avoid an error but it works ok even with error)
 fmax - integer
 strand - 1, 1 or 0 (0 must be passed as string or else will be undef)
 phase - int
 residue_info - string
 locgroup - int (default = 0 so don't pass unless you know its different)
 rank - int (default = 0 so don't pass unless you know its different)
 is_fmin_partial - boolean 't' or 'f' default = 'f'
 is_fmax_partial - boolean 't' or 'f' default = 'f'

=item C<create_ch_featureloc_pub>

 CREATE featureloc_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of featureloc_pub
      or just appending the pub element to featureloc_pub if that is passed
 params
 doc - XML::DOM::Document required
 featureloc_id - OPTIONAL macro featureloc id or XML::DOM featureloc element to create freestanding featureloc_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc (optional featureloc_id) as only params
 uniquename - string required unless pub_id
 type_id - macro id for a pub or XML::DOM cvterm element optional unless creating a new pub 
 type - string from pub type cv                          (i.e. not null value but not part of unique key
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_featureprop>

 CREATE featureprop element
 params
 doc - XML::DOM::Document required
 feature_id - optional macro id for a feature or XML::DOM feature element for standalone featureprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for property type or XML::DOM cvterm element required
 type - string from property type cv 
        Note: will default to making a featureprop from 'property type' cv unless cvname is provided
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_featureprop_pub>

 CREATE featureprop_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of featureprop_pub
      or just appending the pub element to featureprop_pub if that is passed
 params
 doc - XML::DOM::Document required
 featureprop_id  - OPTIONAL macro featureprop id or XML::DOM featureprop element to create freestanding featureprop_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required unless pub_id
 type_id - macro id for a pub or XML::DOM cvterm element optional unless creating a new pub 
 type - string from pub type cv                          (i.e. not null value but not part of unique key
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_genotype>

 CREATE genotype or genotype_id element
 params
 doc - XML::DOM::Document - required
 uniquenname - string required
 description - string optional
 name - string optional
 macro_id - optional string to specify as ID attribute for genotype
 with_id - boolean optional if 1 then genotype_id element is returned

=item C<create_ch_grp>

  CREATE grp element 
 params
 doc - XML::DOM::Document required
 uniquename - string required
 type_id - cvterm macro id or XML::DOM cvterm element required unless type
 type - string for a cvterm Note: will default to using 'SO' cv 
 unless a cvname is provided
 cvname - string optional to specify a cv other than 'SO' for the 
         type_id do not use if providing a cvterm element or macro id to type_id
 name - string optional
 is_analysis - boolean optional default = false
 is_obolete - boolean optional default = false
 macro_id - string optional if provide then add an ID attribute to the top level element of provided value
 with_id - boolean optional if 1 then grp_id element is returned

=item C<create_ch_grpprop>

 CREATE grpprop element
 params
 doc - XML::DOM::Document required
 grp_id - optional macro id for a grp or XML::DOM grp element for standalone grpprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - string from grp property type cv or XML::DOM cvterm element required
        Note: will default to making a from 'grp property type' cv unless cvname is provided
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_grpprop_pub>

 CREATE grpprop_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of grpprop_pub
      or just appending the pub element to grpprop_pub if that is passed
 params
 doc - XML::DOM::Document required
 grpprop_id  - OPTIONAL macro grpprop id or XML::DOM grpprop element to create freestanding grpprop_pub
 pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_grp_relationship>

 CREATE grp_relationship element
 NOTE: this can now create a free standing grp relationship if you pass subject_id or object_id
 params
 doc - XML::DOM::Document required
 object _id - macro id for object grp or XML::DOM grp element
 subject_id - macro id for subject grp or XML::DOM grp element
 NOTE you can pass one or both of the above parameters with the following rules:
    if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other grp
    if both are passed then is_object, is_subject and any parameters to create a grp are ignored
 is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
 is_subject - boolean 't'          this flag indicates if the grp info provided should be 
                                   added in as subject or object grp
 rtype_id - macro id for cvterm or XML::DOM cvterm element
 rtype - string for relationship type note: if relationship name is given will be assigned to 'relationship type' cv
 rank - integer optional with default = 0
 grp_id - macro id for grp or XML::DOM grp element required unless minimal grp bits provided
 uniquename - string required unless grp_id provided
 ftype_id -  macro id for cvterm or XML::DOM cvterm element required unless grp provided
 ftype - string = name of grp type 

=item C<create_ch_grp_relationship_pub>

 CREATE grp_relationship_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of grp_relationship_pub
      or just appending the pub element to grp_relationship_pub if that is passed
 params
 doc - XML::DOM::Document required
 grp_relationship_id -optional grp_relationship XML::DOM element or macro id
 pub_id - XML::DOM pub element - if this is used then pass this and doc (and optional fr_id) as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_grp_relationshipprop>

 CREATE grp_relationshipprop element
 params
 doc - XML::DOM::Document required
 grp_relationship_id - optional macro id for a grp_relationship or XML::DOM grp_relationship element 
                           for standalone grp_relationshipprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for a grp_relationshipprop type or XML::DOM cvterm element required
 type - string from grp_relationshipprop type cv 
        Note: will default to making a grpprop from above cv unless cvname is provided
 cvname - string optional 
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_grp_synonym>

 CREATE grp_synonym element
 params
 doc - XML::DOM::Document required
 grp_id - OPTIONAL macro grp id or XML::DOM grp element to create freestanding grp_synonym
 synonym_id - macro id for synonym or XML::DOM synonym element required unless name and type provided
 name - string required unless synonym_id element provided
 type_id - macro id for cvterm or XML::DOM cvterm element required if name and not type
 type - string = name from the 'synonym type' cv required if name and type_id
 pub_id - macro id for a pub or a XML::DOM pub element required
 pub - a pub uniquename (i.e. FBrf)
 synonym_sgml - string optional but if not provided then synonym_sgml = name
               - do not provide if a synonym_id element is provided
 is_current - optional string = 'f' or 't' default is 't' so don't provide this param
               unless you know you want to change the value
 is_internal - optional string = 't' default is 'f' so don't provide this param
               unless you know you want to change the value

=item C<create_ch_grp_cvterm>

 CREATE grp_cvterm element
 params
 doc - XML::DOM::Document required
 grp_id - OPTIONAL macro grp id or XML::DOM grp element to create freestanding grp_cvterm
 cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
 name - string required unless cvterm_id provided
 cv_id - macro id for cv or XML::DOM cv element required unless cvterm_id provided
 cv - string = cvname
 pub_id - macro id for a pub or XML::DOM pub element required unless pub 
 pub - string pub uniquename Note: cannot create a new pub if uniquename provided
 is_not - optional boolean 't' or 'f' with default = 'f' so don't pass unless you know you want to change

=item C<create_ch_grp_dbxref>

 CREATE grp_dbxref element
 params
 doc - XML::DOM::Document required
 grp_id - OPTIONAL macro grp id or XML::DOM grp element to create freestanding grp_dbxref
 dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
 accession - string required unless dbxref_id provided
 db_id - macro db id or XML::DOM db element required if accession and not db provided
 db - string = db name required if accession and not db_id provided
 version - string optional
 description - string optional
 is_current - string 't' or 'f' boolean default = 't' so don't pass unless
              this should be changed


=item C<create_ch_grp_pub>

 CREATE grp_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of grp_pub
      or just appending the pub element to grp_pub if that is passed 
 params
 doc - XML::DOM::Document required
 grp_id - OPTIONAL macro grp id or XML::DOM grp element to create freestanding grp_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc (optional grp_id)  as only params
 uniquename - string for pub uniquename required unless pub_id
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
           creating a new pub (i.e. not null value but not part of unique key
 type - string from pub type cv same requirement as type_id
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_grp_pubprop>

 CREATE grp_pubprop element
 params
 doc - XML::DOM::Document required
 grp_pub_id - optional macro id for a grp_pub or XML::DOM grp_pub element for standalone grp_pubprop
 value - string - not strictly required
 type_id - macro id for grp_pubprop type or XML::DOM cvterm element required 
 type - string from grp_pubprop type cv
        Note: will default to making a pubprop from above cv unless cvname is provided
 cvname - string optional 
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_grpmember>

 CREATE grpmember element
 params
 doc - XML::DOM::Document required
 grp_id - optional macro id for a grp or XML::DOM grp element for standalone grpmember
 type_id - macro id for grpmember type or XML::DOM cvterm element required 
 type - string from grpmember type cv
        Note: will default to making a grpmember type from above cv unless cvname is provided
 cvname - string optional 
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0
 macro_id - string optional if provide then add an ID attribute to the top level element of provided value
 with_id - boolean optional if 1 then grpmember_id element is returned

=item C<create_ch_grpmember_cvterm>

 CREATE grpmember_cvterm element
 params
 doc - XML::DOM::Document required
 cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
 name - string required unless cvterm_id provided
 cv_id - macro id for cv or XML::DOM cv element required unless cvterm_id provided
 cv - string = cvname
 pub_id - macro id for a pub or XML::DOM pub element required unless pub 
 pub - string pub uniquename Note: cannot create a new pub if uniquename provided
 is_not - optional boolean 't' or 'f' with default = 'f' so don't pass unless you know you want to change

=item C<create_ch_grpmember_pub>

 CREATE grpmember_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of grpmember_pub
      or just appending the pub element to grpmember_pub if that is passed 
 params
 doc - XML::DOM::Document required
 grpmember_id - OPTIONAL macro grpmember id or XML::DOM grpmember element to create freestanding grpmember_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc (optional grpmember_id)  as only params
 uniquename - string for pub uniquename required unless pub_id
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
           creating a new pub (i.e. not null value but not part of unique key
 type - string from pub type cv same requirement as type_id
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_grpmemberprop>

 CREATE grpmemberprop element
 params
 doc - XML::DOM::Document required
 grpmember_id - OPTIONAL macro grpmember id or XML::DOM grpmember element to create freestanding grpmemberprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - string from grpmember property type cv or XML::DOM cvterm element required
        Note: will default to making a from 'grpmember property type' cv unless cvname is provided
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_grpmemberprop_pub>

 CREATE grpmemberprop_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of grpmemberprop_pub
      or just appending the pub element to grpmemberprop_pub if that is passed
 params
 doc - XML::DOM::Document required
 grpmemberprop_id - OPTIONAL macro grpmemberprop id or XML::DOM grpmemberprop element to create freestanding grpmemberprop_pub
 pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_humanhealth>

 CREATE humanhealth element
 params
 doc - XML::DOM::Document required
 uniquename - string required
 name - string optional
 organism_id - organism macro id or XML::DOM organism element required if no genus and species
 genus - string required if no organism
 species - string required if no organism
 dbxref_id - dbxref macro id or XML::DOM dbxref element optional
 accession - optional
 version - optional
 db - optional dbname
 is_obsolete - boolean 't' or 'f' default = 'f' optional
 macro_id - optional string to add as ID attribute value to feature
 with_id - boolean optional if 1 then feature_id element is returned
 no_lookup - boolean optional


=item C<create_ch_humanhealth_cvterm>

 CREATE humanhealth_cvterm element
 params
 doc - XML::DOM::Document required
 humanhealth_id - OPTIONAL macro humanhealth id or XML::DOM humanhealth element to create freestanding humanhealth_cvterm
 cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
 name - string required unless cvterm_id provided
 cv_id - macro id for cv or XML::DOM cv element required unless cvterm_id provided
 cv - string = cvname
 pub_id - macro id for a pub or XML::DOM pub element required unless pub 
 pub - string pub uniquename Note: cannot create a new pub if uniquename provided

=item C<create_ch_humanhealth_cvtermprop>

 CREATE humanhealth_cvtermprop element
 params
 doc - XML::DOM::Document required
 humanhealth_cvterm_id - optional macro id for a humanhealth_cvterm or XML::DOM humanhealth_cvterm element 
                     for standalone humanhealth_cvtermprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for humanhealth_cvtermprop type or XML::DOM cvterm element required
 type -  string from  humanhealth_cvtermprop type cv 
        Note: will default to making a humanhealthprop from 'property type' cv unless cvname is provided
 cvname - string (probably want to pass 'humanhealth_cvtermprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_humanhealth_dbxref>

 CREATE humanhealth_dbxref element
 params
 doc - XML::DOM::Document required
 humanhealth_id - OPTIONAL macro humanhealth id or XML::DOM humanhealth element to create freestanding humanhealth_dbxref
 dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
 accession - string required unless dbxref_id provided
 db_id - macro db id or XML::DOM db element required if accession and not db provided
 db - string = db name required if accession and not db_id provided
 version - string optional
 description - string optional
 macro_id - optional string to add as ID attribute value to humanhealth_dbxref
 is_current - string 't' or 'f' boolean default = 't' so don't pass unless
              this should be changed

=item C <create_ch_humanhealth_dbxrefprop>

 CREATE humanhealth_dbxrefprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for humanhealth_dbxrefprop type or XML::DOM cvterm element required
 type -  string from  property type cv 
 cvname - string (probably want to pass 'humanhealth_dbxrefprop type' if ever used) optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C <create_ch_humanhealth_dbxrefprop_pub>

 CREATE humanhealth_dbxrefprop_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of humanhealthprop_pub
      or just appending the pub element to humanhealthprop_pub if that is passed
 params
 doc - XML::DOM::Document required
 humanhealth_dbxrefprop_id  - OPTIONAL macro humanhealth_dbxrefprop id or XML::DOM humanhealth_dbxrefprop element to create freestanding humanhealth_dbxrefprop_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required unless pub_id
 type_id - macro id for a pub or XML::DOM cvterm element optional unless creating a new pub 
 type - string from pub type cv                          (i.e. not null value but not part of unique key
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_humanhealth_feature>

 CREATE humanhealth_feature element
 params
 doc - XML::DOM::Document required
 organism_id - macro id for organism or XML::DOM organism element
 genus - string
 species - string
 NOTE: you can use the generic parameters in the following cases:
       1.  you are only building either a humanhealth or feature element and not both
       2.  or both humanhealth and feature have the same organism
       otherwise use the prefixed parameters
 WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
 humanhealth_id - macro id for humanhealth or XML::DOM humanhealth element
 hh_uniquename - string humanhealth uniquename
 hh_organism_id - macro id for organism or XML::DOM organism element to link to humanhealth
 hh_genus - string for genus for humanhealth organism
 hh_species - string for species for humanhealth organism
 feature_id - macro id for feature or XML::DOM feature element
 feat_uniquename - string feature uniquename
 feat_organism_id - macro id for organism or XML::DOM organism element to link to feature
 feat_genus - string for genus for feature organism
 feat_species - string for species for feature organism
 feat_type_id - macro id for feature type or XML::DOM cvterm element
 feat_type - string for feature type from SO cv
 pub_id - macro id for pub or XML::DOM pub element
 pub - pub uniquename

=item C <create_ch_humanhealth_featureprop>

 CREATE humanhealth_featureprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for humanhealth_featureprop type or XML::DOM cvterm element required
 type -  string from  humanhealth_featureprop type cv 
        Note: will default to making a featureprop from 'humanhealth_featureprop type' cv unless cvname is provided
 cvname - string (probably want to pass 'humanhealth_featureprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_humanhealth_phenotype>

 CREATE humanhealth_phenotype element
 params
 doc - XML::DOM::Document required
 humanhealth_id - optional macro id for humanhealth or XML::DOM humanhealth element
 phenotype_id - optional macro id for phenotype or XML::DOM phenotype element 
 uniquename - string required if no phenotype_id
 observable_id - macro id for observable or XML::DOM cvterm element optional
 attr_id - macro id for attr or XML::DOM cvterm element optional
 cvalue_id - macro id for cvalue or XML::DOM cvterm element optional
 assay_id - macro id for assay or XML::DOM cvterm element optional
 value - string optional
 macro_id - optional string to specify as ID attribute for phenotype

=item C <create_ch_humanhealth_phenotypeprop>

 CREATE humanhealth_phenotypeprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for humanhealth_phenotypeprop type or XML::DOM cvterm element required
 type -  string from  humanhealth_phenotypeprop type cv 
        Note: will default to making a featureprop from 'humanhealth_phenotypeprop type' cv unless cvname is provided
 cvname - string (probably want to pass 'humanhealth_phenotypeprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

 
=item C<create_ch_humanhealth_pub>

 CREATE humanhealth_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of humanhealth_pub
      or just appending the pub element to humanhealth_pub if that is passed 
 params
 doc - XML::DOM::Document required
 humanhealth_id - OPTIONAL macro humanhealth id or XML::DOM humanhealth element to create freestanding humanhealth_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc (optional humanhealth_id)  as only params
 uniquename - string for pub uniquename required unless pub_id
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
           creating a new pub (i.e. not null value but not part of unique key
 type - string from pub type cv same requirement as type_id
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_humanhealth_pubprop>

 CREATE humanhealth_pubprop element
 params
 doc - XML::DOM::Document required
 humanhealth_pub_id - optional macro id for a humanhealth_pub or XML::DOM humanhealth_pub element for standalone humanhealth_pubprop
 value - string - not strictly required
 type_id - macro id for humanhealth_pubprop type or XML::DOM cvterm element required 
 type - string from humanhealth_pubprop type cv
        Note: will default to making a pubprop from above cv unless cvname is provided
 cvname - string optional 
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_humanhealth_relationship>

 CREATE humanhealth_relationship element
 params
 doc - XML::DOM::Document required
 object _id - macro id for object humanhealth or XML::DOM humanhealth element
 subject_id - macro id for subject humanhealth or XML::DOM humanhealth element
 NOTE you can pass one or both of the above parameters with the following rules:
 if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other humanhealth
 if both are passed then is_object, is_subject and any parameters to create a humanhealth are ignored
 is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
 is_subject - boolean 't'          this flag indicates if the humanhealth info provided should be 
                                   added in as subject or object humanhealth
 rtype_id - macro id for relationship type or XML::DOM cvterm element (Note: currently all is_relationship = '0'
 rtype - string for relationship type note: with this param  type will be assigned to relationship_type cv
 humanhealth_id - macro_id for a humanhealth or XML::DOM humanhealth element required unless minimal humanhealth bits provided
 uniquename - string required unless humanhealth provided
 organism_id - macro id for organism or XML::DOM organism element required unless humanhealth or (genus & species) provided
 genus - string required unless humanhealth or organism provided
 species - string required unless humanhealth or organism provided

=item C<create_ch_humanhealth_relationship_pub>

 CREATE humanhealth_relationship_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of humanhealth_relationship_pub
      or just appending the pub element to humanhealth_relationship_pub if that is passed
 params
 doc - XML::DOM::Document required
 humanhealth_relationship_id -optional humanhealth_relationship XML::DOM element or macro id
 pub_id - XML::DOM pub element - if this is used then pass this and doc (and optional fr_id) as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_humanhealth_synonym>

 CREATE humanhealth_synonym element
 params
 doc - XML::DOM::Document required
 humanhealth_id - OPTIONAL macro humanhealth id or XML::DOM humanhealth element to create freestanding humanhealth_synonym
 synonym_id - macro id for synonym or XML::DOM synonym element required unless name and type provided
 name - string required unless synonym_id element provided
 type_id - macro id for cvterm or XML::DOM cvterm element required if name and not type
 type - string = name from the 'synonym type' cv required if name and type_id
 pub_id - macro id for a pub or a XML::DOM pub element required
 pub - a pub uniquename (i.e. FBrf)
 synonym_sgml - string optional but if not provided then synonym_sgml = name
               - do not provide if a synonym_id element is provided
 is_current - optional string = 'f' or 't' default is 't' so don't provide this param
               unless you know you want to change the value
 is_internal - optional string = 't' default is 'f' so don't provide this param
               unless you know you want to change the value

=item C<create_ch_humanhealthprop>

 CREATE humanhealthprop element
 params
 doc - XML::DOM::Document required
 humanhealth_id - optional macro id for a humanhealth or XML::DOM humanhealth element for standalone humanhealthprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for property type or XML::DOM cvterm element required
 type - string from property type cv 
        Note: will default to making a humanhealthprop from 'property type' cv unless cvname is provided
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_humanhealthprop_pub>

 CREATE humanhealthprop_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of humanhealthprop_pub
      or just appending the pub element to humanhealthprop_pub if that is passed
 params
 doc - XML::DOM::Document required
 humanhealthprop_id  - OPTIONAL macro humanhealthprop id or XML::DOM humanhealthprop element to create freestanding humanhealthprop_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required unless pub_id
 type_id - macro id for a pub or XML::DOM cvterm element optional unless creating a new pub 
 type - string from pub type cv                          (i.e. not null value but not part of unique key
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string



=item C<create_ch_interaction>

 CREATE interaction element
 params
 doc - XML::DOM::Document required
 uniquename - string required
 type_id - cvterm macro id or XML::DOM cvterm element required unless type
 type - string for a cvterm Note: will default to using 'PSI-MI' cv 
 unless a cvname is provided
 cvname - string optional to specify a cv other than 'PSI-MI' for the 
         type_id do not use if providing a cvterm element or macro id to type_id
 description - string optional
 is_obolete - boolean optional default = false
 macro_id - string optional if provide then add an ID attribute to the top level element of provided value
 with_id - boolean optional if 1 then interaction_id element is returned

=item C<create_ch_interactionprop>

 CREATE interactionprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - string from interaction property type cv or XML::DOM cvterm element required
        Note: will default to making a from 'interaction property type' cv unless cvname is provided
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_interactionprop_pub>

 CREATE interactionprop_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of interactionprop_pub
      or just appending the pub element to interactionprop_pub if that is passed
 params
 doc - XML::DOM::Document required
 pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_interaction_cell_line>

 CREATE interaction_cell_line element
 params
 doc - XML::DOM::Document required
 interaction_id - macro id for interaction or XML::DOM interaction element
 int_uniquename - string interaction uniquename
 int_type_id - macro id for interaction type or XML::DOM cvterm element 
 int_type - string for interaction type from interaction type cv
 cell_line_id - macro id for cell_line or XML::DOM cell_line element
 cell_uniquename - string cell_line uniquename
 organism_id - macro id for organism or XML::DOM organism element
 genus - string
 species - string
 NOTE: organism info is only required if you are building a cell_line element
 pub_id -  macro id for pub or XML::DOM pub element
 pub - pub uniquename

=item C<create_ch_interaction_cvterm>

 CREATE interaction_cvterm element
 params
 doc - XML::DOM::Document required
 interaction_id - OPTIONAL macro interaction id or XML::DOM interaction element to create freestanding interaction_cvterm
 cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
 name - string required unless cvterm_id provided
 cv_id - macro id for cv or XML::DOM cv element required unless cvterm_id provided
 cv - string = cvname

=item C<create_ch_interaction_cvtermprop>

 CREATE interaction_cvtermprop element
 params
 doc - XML::DOM::Document required
 interaction_cvterm_id - optional macro id for a interaction_cvterm or XML::DOM interaction_cvterm element 
                     for standalone interaction_cvtermprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for interaction_cvtermprop type or XML::DOM cvterm element required
 type -  string from  interaction_cvtermprop type cv 
        Note: will default to making a interactionprop from 'property type' cv unless cvname is provided
 cvname - string (probably want to pass 'interaction_cvtermprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_interaction_expression>

 CREATE interaction_expression element
 params
 doc - XML::DOM::Document required
 interaction_id - macro id for interaction or XML::DOM interaction element
 int_uniquename - string interaction uniquename
 int_type_id - macro id for interaction type or XML::DOM cvterm element 
 int_type - string for interaction type from interaction type cv
 expression_id - macro id for expression or XML::DOM expression element
 exp_uniquename - string expression uniquename
 pub_id -  macro id for pub or XML::DOM pub element
 pub - pub uniquename

=item C <create_ch_interaction_expressionprop>

 CREATE interaction_expressionprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for interaction_expressionprop type or XML::DOM cvterm element required
 type -  string from  interaction_cvtermprop type cv 
        Note: will default to making a featureprop from 'expression_cvtermprop type' cv unless cvname is provided
 cvname - string (probably want to pass 'interaction_expressionprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_interaction_pub>

 CREATE interaction_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of interaction_pub
      or just appending the pub element to interaction_pub if that is passed
 params
 doc - XML::DOM::Document required
 interaction_id - optional interaction element or macro_id to create freestanding ele
 pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_library>

 CREATE library or library_id element
 params
 doc - XML::DOM::Document
 uniquename - string required
 type_id - macro id for library type or XML::DOM cvterm element required
 type - string for library type from FlyBase miscellaneous CV cv
 organism_id - macro_id for organism or XML::DOM organism element required if no genus and species
 genus - string required if no organism_id
 species - string required if no organism_id
 name - string optional
 macro_id -  optional string to specify as ID attribute for library
 with_id - boolean optional if 1 then db_id element is returned

=item C<create_ch_libraryprop>

 CREATE libraryprop element
 params
 doc - XML::DOM::Document required
 library_id - optional macro id for a library or XML::DOM library element for standalone libraryprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for a library property type or XML::DOM cvterm element required
 type - string from  property type cv
        Note: will default to making a featureprop from 'property type' cv unless cvname is provided
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_library_cvterm>

 CREATE library_cvterm element
 params
 doc - XML::DOM::Document required
 library_id - OPTIONAL macro library id or XML::DOM library element to create freestanding library_cvterm
 cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
 name - string required unless cvterm_id provided
 cv_id - macro id for cv or XML::DOM cv element required unless cvterm_id provided
 cv - string = cvname
 pub_id - macro id for a pub or XML::DOM pub element required unless pub 
 pub - string pub uniquename Note: cannot create a new pub if uniquename provided

=item C<create_ch_library_cvtermprop>

 CREATE library_cvtermprop element
 params
 doc - XML::DOM::Document - require
 library_cvterm_id - optional macro id for a library_cvterm or XML::DOM library_cvterm element 
                        for standalone library_cvtermprop
 value - string
 type_id - macro id for a library_cvtermprop type cvterm or XML::DOM cvterm element required
 type - string from library_cvtermprop type cv
     Note: will default to making a cvtermprop from 'library_cvtermprop type' cv unless cvname is provided
 cvname - string optional but see above for the default 
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C <create_ch_library_dbxref>

 CREATE library_dbxref element
 params
 doc - XML::DOM::Document required
 library_id - OPTIONAL macro library id or XML::DOM library element to create freestanding library_expression 
 dbxref_id - macro dbxref id or XML::DOM dbxref element - required 
 accession - string required unless dbxref_id provided
 db_id - macro db id or XML::DOM db element required if accession and not db provided
 db - string = db name required if accession and not db_id provided
 version - string optional
 description - string optional

=item C <create_ch_library_dbxrefprop>

 CREATE library_dbxrefprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for library_dbxrefprop type or XML::DOM cvterm element required
 type -  string from  property type cv 
 cvname - string (probably want to pass 'library_dbxrefprop type' if ever used) optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C <create_ch_library_expression>

 CREATE library_expression element
 params
 doc - XML::DOM::Document required
 library_id - OPTIONAL macro library id or XML::DOM library element to create freestanding library_expression 
 expression_id - macro expression id or XML::DOM expression element - required unless uniquename provided
 uniquename - string required unless expression_id
 pub_id - macro pub id or XML::DOM pub element - required unless puname provided
 pub - string uniquename for pub (note will have lookup so can't create new pub here)

=item C <create_ch_library_expressionprop>

 CREATE library_expressionprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for library_expressionprop type or XML::DOM cvterm element required
 type -  string from  library_expressionprop type cv 
        Note: will default to making a libraryprop from 'expression_cvtermprop type' cv unless cvname is provided
 cvname - string (probably want to pass 'library_expressionprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_library_feature>

 CREATE library_feature element
 params
 doc - XML::DOM::Document required
 organism_id - macro id for organism or XML::DOM organism element
 genus - string
 species - string
 NOTE: you can use the generic parameters in the following cases:
       1.  you are only building either a library or feature element and not both
       2.  or both library and feature have the same organism
       otherwise use the prefixed parameters
 WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
 library_id - macro id for library or XML::DOM library element
 lib_uniquename - string library uniquename
 lib_organism_id - macro id for organism or XML::DOM organism element to link to library
 lib_genus - string for genus for library organism
 lib_species - string for species for library organism
 lib_type_id - macro id for library type or XML::DOM cvterm element 
 lib_type - string for library type from FlyBase miscellaneous CV cv
 feature_id - macro id for feature or XML::DOM feature element
 feat_uniquename - string feature uniquename
 feat_organism_id - macro id for organism or XML::DOM organism element to link to feature
 feat_genus - string for genus for feature organism
 feat_species - string for species for feature organism
 feat_type_id - macro id for feature type or XML::DOM cvterm element
 feat_type - string for feature type from SO cv

=item C <create_ch_library_featureprop>

 CREATE library_featureprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for library_featureprop type or XML::DOM cvterm element required
 type -  string from  library_featureprop type cv 
        Note: will default to making a featureprop from 'library_featureprop type' cv unless cvname is provided
 cvname - string (probably want to pass 'library_featureprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0


=item C<create_ch_library_grpmember>

 CREATE library_grpmember element
 params
 doc - XML::DOM::Document required
 library_id  optional macro id for a library or XML::DOM library element for standalone library_grpmember
 grpmember_id - macro id for a grpmember or XML::DOM grpmember element for standalone library_grpmember

=item C<create_ch_library_humanhealth>

 CREATE library_humanhealth element
 params
 doc - XML::DOM::Document required
 organism_id - macro id for organism or XML::DOM organism element
 genus - string
 species - string

 NOTE: you can use the generic parameters in the following cases:
       1.  you are only building either a library or humanhealth element and not both
       2.  or both library and humanhealth have the same organism
       otherwise use the prefixed parameters
 WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
 library_id - macro id for library or XML::DOM library element
 lib_uniquename - string library uniquename
 lib_organism_id - macro id for organism or XML::DOM organism element to link to library
 lib_genus - string for genus for library organism
 lib_species - string for species for library organism
 lib_type_id - macro id for library type or XML::DOM cvterm element 
 lib_type - string for library type from FlyBase miscellaneous CV cv
 humanhealth_id - macro id for humanhealth or XML::DOM humanhealth element
 hh_uniquename - string humanhealth uniquename
 hh_organism_id - macro id for organism or XML::DOM organism element to link to humanhealth
 hh_genus - string for genus for humanhealth organism
 hh_species - string for species for humanhealth organism
 pub_id - macro id for pub or XML::DOM pub element
 pub - pub uniquename

=item C <create_ch_library_humanhealthprop>

 CREATE library_humanhealthprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for library_humanhealthprop type or XML::DOM cvterm element required
 type -  string from  library_humanhealthprop type cv 
        Note: will default to making a humanhealthprop from 'library_humanhealthprop type' cv unless cvname is provided
 cvname - string (probably want to pass 'library_humanhealthprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_library_interaction>

 CREATE library_interaction element
 params
 doc - XML::DOM::Document required
 library_id - macro id for library or XML::DOM library element
 lib_uniquename - string library uniquename
 organism_id - macro id for organism or XML::DOM organism element
 genus - string
 species - string
 NOTE: organism info is only required if you are building a library element
 lib_type_id - macro id for library type or XML::DOM cvterm element 
 lib_type - string for library type from FlyBase miscellaneous CV cv
 interaction_id - macro id for interaction or XML::DOM interaction element
 int_uniquename - string interaction uniquename
 int_type_id - macro id for interaction type or XML::DOM cvterm element
 int_type - string for interaction type from psi-mi ontology
 pub_id - macro id for pub or XML::DOM pub element
 pub - pub uniquename

=item C<create_ch_library_pub>

 CREATE library_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of library_pub
      or just appending the pub element to library_pub if that is passed
 params
 doc - XML::DOM::Document required
 library_id - OPTIONAL macro library id or XML::DOM featureloc element to create freestanding library_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this (optional library_id)  and doc as only params
 uniquename - string for pub uniquename required unless pub_id
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
           creating a new pub (i.e. not null value but not part of unique key
 type - string from pub type cv same requirement as type_id
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_library_relationship>

 CREATE library_relationship element
 params
 doc - XML::DOM::Document required
 object _id - macro id for object library or XML::DOM library element
 subject_id - macro id for subject library or XML::DOM library element
 NOTE you can pass one or both of the above parameters with the following rules:
 if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other library
 if both are passed then is_object, is_subject and any parameters to create a library are ignored
 is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
 is_subject - boolean 't'          this flag indicates if the library info provided should be 
                                   added in as subject or object library
 rtype_id - macro id for relationship type or XML::DOM cvterm element (Note: currently all is_relationship = '0'
 rtype - string for relationship type note: with this param  type will be assigned to relationship_type cv
 library_id - macro_id for a library or XML::DOM library element required unless minimal library bits provided
 uniquename - string required unless library provided
 organism_id - macro id for organism or XML::DOM organism element required unless library or (genus & species) provided
 genus - string required unless library or organism provided
 species - string required unless library or organism provided
 ftype_id - macro id for library type or XML::DOM cvterm element required unless library provided
 ftype - string for library type from FlyBase miscellaneous CV

=item C<create_ch_library_relationship_pub>

 CREATE library_relationship_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of library_relationship_pub
      or just appending the pub element to library_relationship_pub if that is passed
 params
 doc - XML::DOM::Document required
 library_relationship_id -optional library_relationship XML::DOM element or macro id
 pub_id - XML::DOM pub element - if this is used then pass this and doc (and optional fr_id) as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_library_strain>

 CREATE library_strain element
 params
 doc - XML::DOM::Document required
 organism_id - macro id for organism or XML::DOM organism element
 genus - string
 species - string

 NOTE: you can use the generic parameters in the following cases:
       1.  you are only building either a library or strain element and not both
       2.  or both library and strain have the same organism
       otherwise use the prefixed parameters
 WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
 library_id - macro id for library or XML::DOM library element
 lib_uniquename - string library uniquename
 lib_organism_id - macro id for organism or XML::DOM organism element to link to library
 lib_genus - string for genus for library organism
 lib_species - string for species for library organism
 lib_type_id - macro id for library type or XML::DOM cvterm element 
 lib_type - string for library type from FlyBase miscellaneous CV cv
 strain_id - macro id for strain or XML::DOM strain element
 str_uniquename - string strain uniquename
 str_organism_id - macro id for organism or XML::DOM organism element to link to strain
 str_genus - string for genus for strain organism
 str_species - string for species for strain organism
 pub_id - macro id for pub or XML::DOM pub element
 pub - pub uniquename

=item C <create_ch_library_strainprop>

 CREATE library_strainprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for library_strainprop type or XML::DOM cvterm element required
 type -  string from  library_strainprop type cv 
        Note: will default to making a strainprop from 'library_strainprop type' cv unless cvname is provided
 cvname - string (probably want to pass 'library_strainprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0


=item C<create_ch_library_synonym>

 CREATE library_synonym element
 params
 doc - XML::DOM::Document required
 library_id - OPTIONAL macro library id or XML::DOM library element to create freestanding library_synonym
 synonym_id - macro id for synonym or XML::DOM synonym element required unless name and type provided
 name - string required unless synonym_id element provided
 type_id - macro id for cvterm or XML::DOM cvterm element required if name and not type
 type - string = name from the 'synonym type' cv required if name and type_id
 pub_id - macro id for a pub or a XML::DOM pub element required
 pub - a pub uniquename (i.e. FBrf)
 synonym_sgml - string optional but if not provided then synonym_sgml = name
               - do not provide if a synonym element is provided
 is_current - optional string = 'f' or 't' default is 't' so don't provide this param
               unless you know you want to change the value
 is_internal - optional string = 't' default is 'f' so don't provide this param
               unless you know you want to change the value

=item C<create_ch_organism>

 CREATE organism or organism_id element
 params
 doc - XML::DOM::Document
 genus - string required
 species - string required
 abbreviation - string optional
 common_name - string optional
 comment - string optional
 macro_id - optional string to add as ID attribute value to organism
 with_id - boolean optional if 1 then organism_id element is returned

=item C<create_ch_organism_cvterm>

 CREATE organism_cvterm 
 params
 doc - XML::DOM::Document - required
 organism_id - OPTIONAL macro organism id or XML::DOM organism element to create freestanding organism_cvterm
 cvterm_id - macro id for a cvterm or XML::DOM cvterm element - required unless name and cv params
 name - string name for cvterm required unless cvterm_id
 cv_id - macro id for cv or XML::DOM cv element required unless cvterm_id or cv provided
 cv - string = cvname required unless cvterm_id or cv_id provided
 rank - integer with default of 0 - this plus type_id used for ordering items 
 pub_id - macro id for pub or XML::DOM pub element required unless pub provided
 pub - string pub uniquename

=item C<create_ch_organism_cvtermprop>

 CREATE organism_cvtermprop
 params
 doc - XML::DOM::Document - require
 organism_cvterm_id - optional macro id for a organism_cvterm or XML::DOM organism_cvterm element 
                        for standalone organism_cvtermprop
 value - string
 type_id - macro id for a organism_cvtermprop type cvterm or XML::DOM cvterm element required
 type - string from organism_cvtermprop type cv
     Note: will default to making a cvtermprop from 'organism_cvtermprop type' cv unless cvname is provided
 cvname - string optional but see above for the default 
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_organism_dbxref>

 CREATE organism_dbxref element
 params
 doc - XML::DOM::Document required
 organism_id - OPTIONAL macro feature id or XML::DOM feature element to create freestanding feature_dbxref
 dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
 accession - string required unless dbxref_id provided
 db_id - macro db id or XML::DOM db element required if accession and not db provided
 db - string = db name required if accession and not db_id provided
 version - string optional
 description - string optional
 is_current - string 't' or 'f' boolean default = 't' so don't pass unless
              this should be changed

=item C<create_ch_organism_grpmember>

 CREATE organism_grpmember element
 params
 doc - XML::DOM::Document required
 organism_id  macro id for a organism or XML::DOM organism element for standalone organism_grpmember
 grpmember_id - macro id for a grpmember or XML::DOM grpmember element for standalone organism_grpmember

=item C<create_ch_organism_library>

 CREATE organism_library element
 params
 doc - XML::DOM::Document required
 organism_id  macro id for a organism or XML::DOM organism element for standalone organism_library
 library_id  optional macro id for a library or XML::DOM library element for standalone organism_library

=item C<create_ch_organism_pub>

 CREATE organism_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of organism_pub
      or just appending the pub element to organism_pub if that is passed 
 params
 doc - XML::DOM::Document required
 organism_id - OPTIONAL macro organism id or XML::DOM organism element to create freestanding organism_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc (optional organism_id)  as only params
 uniquename - string for pub uniquename required unless pub_id
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
           creating a new pub (i.e. not null value but not part of unique key
 type - string from pub type cv same requirement as type_id
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_organismprop>

 CREATE organismprop element
 params
 doc - XML::DOM::Document required
 organism_id - optional macro id for a organism or XML::DOM organism element for standalone organismprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for property type or XML::DOM cvterm element required
 type - string from organismprop type cv 
        Note: will default to making a featureprop from 'organismprop type' cv unless cvname is provided
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_organismprop_pub>

 CREATE organismprop_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of organismprop_pub
      or just appending the pub element to organismprop_pub if that is passed
 params
 doc - XML::DOM::Document required
 organismprop_id -optional organismprop XML::DOM element or macro id
 pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_phendesc>

 CREATE phendesc element
 params
 doc - XML::DOM::Document required
 genotype_id - macro id for genotype or XML::DOM genotype element
 genotype - string genotype uniquename 
 environment_id - macro id for environment or XML::DOM environment element
 environment - string environment uniquename
 description - string optional but if creating a new phendesc this can't be null
 type_id - macro id for phendesc type or XML::DOM cvterm element
 type - string for cvterm name from phendesc type CV
 pub_id - macro id for pub or XML::DOM pub element
 pub - string pub uniquename

=item C<create_ch_phenotype>

 CREATE phenotype or phenotype_id element
 params
 doc - XML::DOM::Document required
 uniquename - string required
 observable_id - macro id for observable or XML::DOM cvterm element optional
 attr_id - macro id for attr or XML::DOM cvterm element optional
 cvalue_id - macro id for cvalue or XML::DOM cvterm element optional
 assay_id - macro id for assay or XML::DOM cvterm element optional
 value - string optional
 macro_id - optional string to specify as ID attribute for genotype
 with_id - boolean optional if 1 then genotype_id element is returned

=item C<create_ch_phenotype_comparison>

 CREATE phenotype_comparison element
 params
 doc - XML::DOM::Document required
 organism_id - macro id for an organism or XML::DOM organism element required unless genus and species
 genus - string required if no organism_id
 species - string required if no organism_id
 genotype1_id - macro id for a genotype or XML::DOM genotype element required unless genotype1
 genotype1 - string genotype uniquename required unless genotype1_id
 environment1_id - macro id for a environment or XML::DOM environment element required unless environment1
 environment1 - string environment uniquename required unless environment1_id
 genotype2_id - macro id for a genotype or XML::DOM genotype element required unless genotype2
 genotype2 - string genotype uniquename required unless genotype2_id
 environment2_id - macro id for a environment or XML::DOM environment element required unless environment2
 environment2 - string environment uniquename required unless environment2_id
 phenotype1_id - macro id for phenotype or XML::DOM phenotype element required unless phenotype1
 phenotype1 - string phenotype uniquename required unless phenotype1_id
 phenotype2_id - macro id for phenotype or XML::DOM phenotype element optional
 phenotype2 - string phenotype uniquename optional
 pub_id macro id for a pub or a XML::DOM pub element required unless pub
 pub - a pub uniquename (i.e. FBrf) required unless pub_id

 Alias: create_ch_ph_comp

=item C<create_ch_phenotype_comparison_cvterm>

 CREATE phenotype_comparison_cvterm element
 params
 doc - XML::DOM::Document optional - required
 phenotype_comparison_id - optional macro id for phenotype_comparison or phenotype_comparison XML::DOM element
 NOTE: to make a standalone element
 cvterm_id -  macro id for cvterm of XML::DOM cvterm element
 name - cvterm name
 cv_id - macro id for a CV or XML::DOM cv element
 cv - name of a cv
 is_obsolete - optional param for cvterm
 NOTE: you need to pass cvterm bits if attaching to existing phenotype element or 
       creating a freestanding phenotype_cvterm
 rank - optional with default = 0 so only pass if you want a different rank

 Alias: create_ch_ph_comp_cvt

=item C<create_ch_phenotype_cvterm>

 CREATE phenotype_cvterm element
 params
 doc - XML::DOM::Document optional - required
 phenotype_id - macro id for phenotype or XML::DOM phenotype element
 uniquename - phenotype uniquename
 NOTE: you need to pass phenotype bits if attaching to existing cvterm element or 
       creating a freestanding phenotype_cvterm
 cvterm_id -  macro id for cvterm of XML::DOM cvterm element
 name - cvterm name
 cv_id - macro id for a CV or XML::DOM cv element
 cv - name of a cv
 is_obsolete - optional param for cvterm
 NOTE: you need to pass cvterm bits if attaching to existing phenotype element or 
       creating a freestanding phenotype_cvterm
 rank - optional with default = 0 so only pass if you want a different rank

 Alias: create_ch_ph_cvt

=item C<create_ch_phenstatement>

 CREATE phenstatement element
 params
 doc - XML::DOM::Document required
 genotype_id - macro id for a genotype or XML::DOM genotype element required unless genotype
 genotype - string genotype uniquename required unless genotype_id
 environment_id - macro id for a environment or XML::DOM environment element required unless environment
 environment - string environment uniquename required unless environment_id
 phenotype_id - macro id for phenotype or XML::DOM phenotype element required unless phenotype
 phenotype - string phenotype uniquename required unless phenotype_id
 type_id - macro id for a phenstatement type or XML::DOM cvterm element
 pub_id - macro id for a pub or a XML::DOM pub element required unless pub
 pub - a pub uniquename (i.e. FBrf) required unless pub_id

=item C<create_ch_pub>

 CREATE pub or pub_id element
 params
 doc - XML::DOM::Document required
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
        creating a new pub (i.e. not null value but not part of unique key
 type - string from pub type cv
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string
 macro_id - optional string to add as ID attribute value to pub
 with_id - boolean optional if 1 then pub_id element is returned
 no_lookup - boolean optional

=item C<create_ch_pubauthor>

 CREATE pubauthor element
 params
 doc - XML::DOM::Document required
 pub_id -  macro pub id or XML::DOM pub element optional to create a freestanding pubauthor element
 pub - pub uniquename optional but required if making a freestanding element unless pub_id 
 rank - positive integer required 
 surname - string - required if creating a pubauthor element but optional for other operations
 editor - boolean 't' or 'f' default = 'f' so don't pass unless you want to change
 givennames - string optional
 suffix - string optional  

=item C<create_ch_pubprop>

 CREATE pubprop element
 params
 doc - XML::DOM::Document required
 pub_id - optional macro id for a pub or XML::DOM pub element for standalone pubprop
 value - string - not strictly required and in some cases this value is null in chado
 type_id - macro id for a pubprop type or XML::DOM cvterm element required
 type - string from pubprop type cv 
        Note: will default to making a pubprop from above cv unless cvname is provided
 cvname - string optional 
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_pub_dbxref>

 CREATE pub_dbxref element
 params
 doc - XML::DOM::Document required
 pub_id - macro pub id or XML::DOM pub element optional to create freestanding pub_dbxref
 dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
 accession - string required unless dbxref_id provided
 db_id - macro db id or XML::DOM db element required unless dbxref_id provided
 db - string name of db
 is_current - string 't' or 'f' boolean default = 't' so don't pass unless

=item C<create_ch_pub_relationship>

 CREATE pub_relationship element
NOTE: this can now create a free standing pub relationship if you pass subject_id or object_id
 params
 doc - XML::DOM::Document required
 object _id - macro id for object feature or XML::DOM feature element
 subject_id - macro id for subject feature or XML::DOM feature element
 is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
 is_subject - boolean 't'          this flag indicates if the pub info provided should be 
                                   added in as subject or object pub
 rtype_id -  macro for relationship type or XML::DOM cvterm element 
 rtype - string for relationship type note: if relationship name is given will be assigned to relationship_type cv
 pub_id - macro id for a pub or XML::DOM pub element required unless uniquename provided
 uniquename - uniquename of the pub - required unless pub element provided
 type_id - macro id for a pub type or XML::DOM cvterm element for pub type
 type - string specifying pub type

 Alias: create_ch_pr

=item C<create_ch_strain>

 CREATE strain element
 params
 doc - XML::DOM::Document required
 uniquename - string required
 name - string optional
 organism_id - organism macro id or XML::DOM organism element required if no genus and species
 genus - string required if no organism
 species - string required if no organism
 dbxref_id - dbxref macro id or XML::DOM dbxref element optional
 accession - optional
 version - optional
 db - optional dbname
 is_obsolete - boolean 't' or 'f' default = 'f' optional
 macro_id - optional string to add as ID attribute value to feature
 with_id - boolean optional if 1 then feature_id element is returned
 no_lookup - boolean optional

=item C<create_ch_strain_cvterm>

 CREATE strain_cvterm element
 params
 doc - XML::DOM::Document required
 strain_id - OPTIONAL macro strain id or XML::DOM strain element to create freestanding strain_cvterm
 cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
 name - string required unless cvterm_id provided
 cv_id - macro id for cv or XML::DOM cv element required unless cvterm_id provided
 cv - string = cvname
 pub_id - macro id for a pub or XML::DOM pub element required unless pub 
 pub - string pub uniquename Note: cannot create a new pub if uniquename provided

=item C<create_ch_strain_cvtermprop>

 CREATE strain_cvtermprop element
 params
 doc - XML::DOM::Document required
 strain_cvterm_id - optional macro id for a strain_cvterm or XML::DOM strain_cvterm element 
                     for standalone strain_cvtermprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for strain_cvtermprop type or XML::DOM cvterm element required
 type -  string from  strain_cvtermprop type cv 
        Note: will default to making a strainprop from 'property type' cv unless cvname is provided
 cvname - string (probably want to pass 'strain_cvtermprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_strain_dbxref>

 CREATE strain_dbxref element
 params
 doc - XML::DOM::Document required
 strain_id - OPTIONAL macro strain id or XML::DOM strain element to create freestanding strain_dbxref
 dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
 accession - string required unless dbxref_id provided
 db_id - macro db id or XML::DOM db element required if accession and not db provided
 db - string = db name required if accession and not db_id provided
 version - string optional
 description - string optional
 is_current - string 't' or 'f' boolean default = 't' so don't pass unless
              this should be changed

=item C<create_ch_strain_feature>

 CREATE strain_feature element
 params
 doc - XML::DOM::Document required
 organism_id - macro id for organism or XML::DOM organism element
 genus - string
 species - string
 NOTE: you can use the generic parameters in the following cases:
       1.  you are only building either a strain or feature element and not both
       2.  or both strain and feature have the same organism
       otherwise use the prefixed parameters
 WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
 strain_id - macro id for strain or XML::DOM strain element
 str_uniquename - string strain uniquename
 str_organism_id - macro id for organism or XML::DOM organism element to link to strain
 str_genus - string for genus for strain organism
 str_species - string for species for strain organism
 str_type_id - macro id for strain type or XML::DOM cvterm element 
 str_type - string for strain type from strain type cv
 feature_id - macro id for feature or XML::DOM feature element
 feat_uniquename - string feature uniquename
 feat_organism_id - macro id for organism or XML::DOM organism element to link to feature
 feat_genus - string for genus for feature organism
 feat_species - string for species for feature organism
 feat_type_id - macro id for feature type or XML::DOM cvterm element
 feat_type - string for feature type from SO cv
 pub_id - macro id for pub or XML::DOM pub element
 pub - pub uniquename

=item C <create_ch_strain_featureprop>

 CREATE strain_featureprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for strain_featureprop type or XML::DOM cvterm element required
 type -  string from  strain_featureprop type cv 
        Note: will default to making a featureprop from 'strain_featureprop type' cv unless cvname is provided
 cvname - string (probably want to pass 'strain_featureprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_strain_phenotype>

 CREATE strain_phenotype element
 params
 doc - XML::DOM::Document required
 strain_id - optional macro id for strain or XML::DOM strain element
 phenotype_id - optional macro id for phenotype or XML::DOM phenotype element 
 uniquename - string required if no phenotype_id
 observable_id - macro id for observable or XML::DOM cvterm element optional
 attr_id - macro id for attr or XML::DOM cvterm element optional
 cvalue_id - macro id for cvalue or XML::DOM cvterm element optional
 assay_id - macro id for assay or XML::DOM cvterm element optional
 value - string optional
 macro_id - optional string to specify as ID attribute for phenotype

=item C <create_ch_strain_phenotypeprop>

 CREATE strain_phenotypeprop element
 params
 doc - XML::DOM::Document required
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for strain_phenotypeprop type or XML::DOM cvterm element required
 type -  string from  strain_phenotypeprop type cv 
        Note: will default to making a featureprop from 'strain_phenotypeprop type' cv unless cvname is provided
 cvname - string (probably want to pass 'strain_phenotypeprop type') optional
          but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

 
=item C<create_ch_strain_pub>

 CREATE strain_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of strain_pub
      or just appending the pub element to strain_pub if that is passed 
 params
 doc - XML::DOM::Document required
 strain_id - OPTIONAL macro strain id or XML::DOM strain element to create freestanding strain_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc (optional strain_id)  as only params
 uniquename - string for pub uniquename required unless pub_id
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
           creating a new pub (i.e. not null value but not part of unique key
 type - string from pub type cv same requirement as type_id
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_strain_relationship>

 CREATE strain_relationship element
 params
 doc - XML::DOM::Document required
 object _id - macro id for object strain or XML::DOM strain element
 subject_id - macro id for subject strain or XML::DOM strain element
 NOTE you can pass one or both of the above parameters with the following rules:
 if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other strain
 if both are passed then is_object, is_subject and any parameters to create a strain are ignored
 is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
 is_subject - boolean 't'          this flag indicates if the strain info provided should be 
                                   added in as subject or object strain
 rtype_id - macro id for relationship type or XML::DOM cvterm element (Note: currently all is_relationship = '0'
 rtype - string for relationship type note: with this param  type will be assigned to relationship_type cv
 strain_id - macro_id for a strain or XML::DOM strain element required unless minimal strain bits provided
 uniquename - string required unless strain provided
 organism_id - macro id for organism or XML::DOM organism element required unless strain or (genus & species) provided
 genus - string required unless strain or organism provided
 species - string required unless strain or organism provided

=item C<create_ch_strain_relationship_pub>

 CREATE strain_relationship_pub element
 Note that this is just calling create_ch_pub 
      and adding returned pub_id element as a child of strain_relationship_pub
      or just appending the pub element to strain_relationship_pub if that is passed
 params
 doc - XML::DOM::Document required
 strain_relationship_id -optional strain_relationship XML::DOM element or macro id
 pub_id - XML::DOM pub element - if this is used then pass this and doc (and optional fr_id) as only params
 uniquename - string required
 type_id - macro id for pub type or XML::DOM cvterm element optional unless 
 type -  string from pub type cv 
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_strain_synonym>

 CREATE strain_synonym element
 params
 doc - XML::DOM::Document required
 strain_id - OPTIONAL macro strain id or XML::DOM strain element to create freestanding strain_synonym
 synonym_id - macro id for synonym or XML::DOM synonym element required unless name and type provided
 name - string required unless synonym_id element provided
 type_id - macro id for cvterm or XML::DOM cvterm element required if name and not type
 type - string = name from the 'synonym type' cv required if name and type_id
 pub_id - macro id for a pub or a XML::DOM pub element required
 pub - a pub uniquename (i.e. FBrf)
 synonym_sgml - string optional but if not provided then synonym_sgml = name
               - do not provide if a synonym_id element is provided
 is_current - optional string = 'f' or 't' default is 't' so don't provide this param
               unless you know you want to change the value
 is_internal - optional string = 't' default is 'f' so don't provide this param
               unless you know you want to change the value

=item C<create_ch_strainprop>

 CREATE strainprop element
 params
 doc - XML::DOM::Document required
 strain_id - optional macro id for a strain or XML::DOM strain element for standalone strainprop
 value - string - not strictly required but if you don't provide this then not much point
 type_id - macro id for property type or XML::DOM cvterm element required
 type - string from property type cv 
        Note: will default to making a strainprop from 'property type' cv unless cvname is provided
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

=item C<create_ch_strainprop_pub>

 CREATE strainprop_pub element
 Note that this is just calling create_ch_pub setting with_id = 1 
      and adding returned pub_id element as a child of strainprop_pub
      or just appending the pub element to strainprop_pub if that is passed
 params
 doc - XML::DOM::Document required
 strainprop_id  - OPTIONAL macro strainprop id or XML::DOM strainprop element to create freestanding strainprop_pub
 pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc as only params
 uniquename - string required unless pub_id
 type_id - macro id for a pub or XML::DOM cvterm element optional unless creating a new pub 
 type - string from pub type cv                          (i.e. not null value but not part of unique key
 title - optional string
 volumetitle - optional string
 volume - optional string
 series_name - optional string
 issue - optional string
 pyear - optional string
 pages - optional string
 miniref - optional string
 is_obsolete - optional string 't' boolean value with default = 'f'
 publisher - optional string

=item C<create_ch_synonym>

 CREATE synonym or synonym_id element
 params
 doc - XML::DOM::Document required
 name - string required
 synonym_sgml - string optional but if not provided then synonym_sgml = name
 type_id - macro id for synonym type or XML::DOM cvterm element
 type - string = name from the 'synonym type' cv
 macro_id - optional string to add as ID attribute value to synonym
 with_id - boolean optional if 1 then synonym_id element is returned

=item C<create_ch_prop>

 GENERIC METHOD FOR CREATING ANY TYPE OF PROP element
 params
 doc - XML::DOM::Document required
 parentname - string that is the name of the element that you want to attach the prop element to
              eg. pass 'feature' to make 'featureprop' element
 parent_id - macro id for the parent table or XML::DOM table_id element
             NOTE: name of this parameter should match table_id for 
             type of prop eg. feature_id for a featureprop or pub_id for pubprop
 value - string - not strictly required but if you don't provide this then not much point actually
 type_id - macro id for a property type or XML::DOM cvterm element required
 type - string from a property type cv 
           Note: will default to making a type of from 'tablenameprop type' cv unless cvname is provided
  WARNING: as property type cv names are not consistent SAFEST to provide cvname
 cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
 rank - integer with a default of zero so don't use unless you want a rank other than 0

 NOTE: this method is now called by all the create_ch_xxxxprop methods which are really just
       wrappers that provide the most likely desired cvname for the property type unless you
       provide one

=item C< _add_orgid_param>

 internal helper method to swap parameters from genus species to organism_id

=item C<_create_simple_element>

 internal helper method that is called by functions to build simple table elements
 i.e. those that do not reference another table
 eg. contact, expression, genotype, organism

=item C<_build_element>

 internal helper method to build up elements

=back

=head1 AUTHOR

Andy Schroeder - andy@morgan.harvard.edu

=head1 SEE ALSO

PrettyPrintDom,  XML::DOM

=cut

# CREATE analysis or analysis_id element
# params 
# doc - XML::DOM::Document optional - required
# program - required string
# programversion - required string (usually a number) NOTE: can add default 1.0?
# sourcename - optional string NOTE: this is part of the unique key and while no constraint usually shouldn't be null
# name - optional string
# description - optional string
# algorithm - optional string
# sourceversion - optional string
# sourceuri -optional string
# timeexecuted - optional string value like '1999-01-08 04:05:06' default will be whenever data is added
#                NOTE: not sure on xort-postgres interaction regarding invalid timestamp formats
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value
# with_id - optional if true will create analysis_id at top level
sub create_ch_analysis {
    my %params = @_;
    print "WARNING -- While 'sourcename' is not required it is usually useful to provide this and you haven't\n"
      unless $params{sourcename};
    $params{required} = ['program','programversion'];
    $params{elname} = 'analysis';
    my $ael = _create_simple_element(%params);
    return $ael;
}

# CREATE analysisfeature element
# params 
# doc - XML::DOM::Document optional - required
#
# Here are parameters to make a feature element must either pass a feature_id or the other necessary bits
# feature_id - macro feature id or XML::DOM feature element
# uniquename - string
# organism_id - macro organism id or XML::DOM organism element
# genus - string
# species - string
# type_id -  macro id for feature type or XML::DOM cvterm element for type
# type - string valid SO feature type
#
# Here are parameters to make an analysis element must either pass analysis_id or required bits
# analysis_id - macro analysis id or XML::DOM analysis element
# program - string
# programversion - string
# sourcename - string
#
# Here are the optional bits that can be added to the analysisfeature
# rawscore - number (double) 
# normscore - number (double)  
# significance - number (double)
# identity - number (double)    
sub create_ch_analysisfeature {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $af_el = $ldoc->createElement('analysisfeature');

    # this first section identifies if a feature needs to be dealt with
    #create a feature element if params are provided
    if ($params{uniquename}) {
      print "You don't have all the parameters required to make a feature, NO GO!\n" and return
	unless ($params{organism_id} or ($params{genus} and $params{species})) and ($params{type_id} or $params{type});

      my @f_ok = qw(doc uniquename organism_id genus species type_id type); # valid feature parameters
      my %fparams;
      # populate the feature parameter hash 
      foreach my $p (keys %params) {
	if (grep $_ eq $p, @f_ok) {
	  $fparams{$p} = $params{$p};
	  delete $params{$p} unless $p eq 'doc';
	}
      }
      $params{feature_id} = create_ch_feature(%fparams);
    }
      
    # or if a macro id or existing element have been provided
    if ($params{feature_id}) {
      $af_el->appendChild(_build_element($ldoc,'feature_id',$params{feature_id}));
      delete $params{feature_id};
    }


    # and here we are dealing with analysis info if provided
    if ($params{program} or $params{programversion} or $params{sourcename}) {
      print "WARNING -- You are trying to make an analysis without providing both program and programversion, NO GO!\n"
	unless ($params{program} and $params{programversion});
      
      my @a_ok = qw(doc program programversion sourcename); # valid feature parameters
      my %aparams;
      # populate the feature parameter hash 
      foreach my $p (keys %params) {
	if (grep $_ eq $p, @a_ok) {
	  $aparams{$p} = $params{$p};
	  delete $params{$p} unless $p eq 'doc';
	}
      }
      $params{analysis_id} = create_ch_analysis(%aparams);
    }

    if ($params{analysis_id}) {
      $af_el->appendChild(_build_element($ldoc,'analysis_id',$params{analysis_id}));
      delete $params{analysis_id};
    }


    foreach my $e (keys %params) {
      next if ($e eq 'doc');
      print "WARNING -- $e should be a valid double number and it's not\n" 
	if $params{$e} !~  /[+-]?(\d+\.\d+|\d+\.|\.\d+)/;  # floating point, no exponent
      $af_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }    

    return $af_el;
}

# CREATE analysisprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from realtionship property type cv or XML::DOM cvterm element required 
#        Note: will default to making a featureprop from 'analysis property type' cv unless cvname is provided
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_analysisprop {
    my %params = @_;
    $params{parentname} = 'analysis';
    unless ($params{type_id}) {
	$params{cvname} = 'analysis property type' unless $params{cvname};
    }
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE analysisgrp element
# params 
# doc - XML::DOM::Document optional - required
#
# Here are parameters to make a grp element must either pass a grp_id or the other necessary bits
# grp_id - macro grp id or XML::DOM grp element
# uniquename - string
# type_id -  macro id for grp type or XML::DOM cvterm element for type
# type - string valid SO grp type
#
# Here are parameters to make an analysis element must either pass analysis_id or required bits
# analysis_id - macro analysis id or XML::DOM analysis element
# program - string
# programversion - string
# sourcename - string
#
# Here are the optional bits that can be added to the analysisgrp
# rawscore - number (double) 
# normscore - number (double)  
# significance - number (double)
# identity - number (double)    
sub create_ch_analysisgrp {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $af_el = $ldoc->createElement('analysisgrp');

    # this first section identifies if a grp needs to be dealt with
    #create a feature element if params are provided
    if ($params{uniquename}) {
      print "You don't have all the parameters required to make a grp, NO GO!\n" and return
	unless ($params{type_id} or $params{type});

      my @f_ok = qw(doc uniquename type_id type); # valid grp parameters
      my %fparams;
      # populate the grp parameter hash 
      foreach my $p (keys %params) {
	if (grep $_ eq $p, @f_ok) {
	  $fparams{$p} = $params{$p};
	  delete $params{$p} unless $p eq 'doc';
	}
      }
      $params{grp_id} = create_ch_grp(%fparams);
    }
      
    # or if a macro id or existing element have been provided
    if ($params{grp_id}) {
      $af_el->appendChild(_build_element($ldoc,'grp_id',$params{grp_id}));
      delete $params{grp_id};
    }

   # and here we are dealing with analysis info if provided
    if ($params{program} or $params{programversion} or $params{sourcename}) {
      print "WARNING -- You are trying to make an analysis without providing both program and programversion, NO GO!\n"
	unless ($params{program} and $params{programversion});
      
      my @a_ok = qw(doc program programversion sourcename); # valid analysis parameters
      my %aparams;
      # populate the analysis parameter hash 
      foreach my $p (keys %params) {
	if (grep $_ eq $p, @a_ok) {
	  $aparams{$p} = $params{$p};
	  delete $params{$p} unless $p eq 'doc';
	}
      }
      $params{analysis_id} = create_ch_analysis(%aparams);
    }

    if ($params{analysis_id}) {
      $af_el->appendChild(_build_element($ldoc,'analysis_id',$params{analysis_id}));
      delete $params{analysis_id};
    }


    foreach my $e (keys %params) {
      next if ($e eq 'doc');
      print "WARNING -- $e should be a valid double number and it's not\n" 
	if $params{$e} !~  /[+-]?(\d+\.\d+|\d+\.|\.\d+)/;  # floating point, no exponent
      $af_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }    

    return $af_el;
}

# CREATE analysisgrpmember element
# params 
# doc - XML::DOM::Document optional - required
#
# Here are parameters to make a grpmember element must either pass a grpmember_id or the other necessary bits
# grpmember_id - macro grpmember id or XML::DOM grpmember element
#
# Here are parameters to make an analysis element must either pass analysis_id or required bits
# analysis_id - macro analysis id or XML::DOM analysis element

sub create_ch_analysisgrpmember {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $af_el = $ldoc->createElement('analysisgrpmember');
     
    # or if a macro id or existing element have been provided
    if ($params{grpmember_id}) {
      $af_el->appendChild(_build_element($ldoc,'grpmember_id',$params{grpmember_id}));
      delete $params{grpmember_id};
    }

    if ($params{analysis_id}) {
      $af_el->appendChild(_build_element($ldoc,'analysis_id',$params{analysis_id}));
      delete $params{analysis_id};
    }

    return $af_el;
}

# CREATE cell_line element
# params
# doc - XML::DOM::Document required
# uniquename - string required
# organism_id - organism macro id or XML::DOM organism element required if no genus and species
# genus - string required if no organism
# species - string required if no organism
# name - string optional
# with_id - boolean optional if 1 then cell_line_id element is returned
# no_lookup - boolean option if 1 then default op="lookup" attribute will not be added to element
sub create_ch_cell_line {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fid_el = $ldoc->createElement('cell_line_id') if $params{with_id};

    ## cell_line element (will be returned)
    my $f_el = $ldoc->createElement('cell_line');
    $f_el->setAttribute('id',$params{macro_id}) if $params{macro_id};    

    # add an op="lookup" attribute unless no_lookup is specified
    unless ($params{no_lookup}) {
	$f_el->setAttribute('op','lookup');
    }
	
    #create organism_id element if genus and species are provided
    unless ($params{organism_id}) {
	_add_orgid_param(\%params);
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id' || $e eq 'no_lookup');
	$f_el->appendChild(_build_element($ldoc, $e,$params{$e}));
    }

    if ($fid_el) {
	$fid_el->appendChild($f_el);
	return $fid_el;
    }
    return $f_el;
}

# CREATE cell_lineprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from  cell_lineprop type cv or XML::DOM cvterm element required 
#        Note: will default to making a featureprop from 'cell_lineprop type' cv unless 
#              cvname is provided
# cvname - string optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_cell_lineprop {
    my %params = @_;
    $params{parentname} = 'cell_line';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE cell_lineprop_pub element
# Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of cell_lineprop_pub
#      or just appending the pub element to featureprop_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_cell_lineprop_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('cell_lineprop_pub');

    if ($params{cell_lineprop_id}) {
	$fp_el->appendChild(_build_element($ldoc,'cell_lineprop_id',$params{cell_lineprop_id}));
	delete $params{cell_lineprop_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE cell_line_cvterm element
# params
# doc - XML::DOM::Document required
# cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
# name - string required unless cvterm_id provided Note: a cvterm has a lookup by default cannot make a new cvterm 
#                                                        with this method
# cv_id - macro id for cv or XML::DOM cv element required if name and not cv
# cv - string for name of cv required if name and not cv_id
# pub_id - macro id for pub or XML::DOM pub element required unless pub
# pub - string = pub uniquename Note: as pub has lookup option by default can't make a new pub using this param
# rank - int optional default = 0
sub create_ch_cell_line_cvterm {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fct_el = $ldoc->createElement('cell_line_cvterm');

    #create a cvterm element if necessary
    unless ($params{cvterm_id}) {
      print "WARNING -- you haven't provided params to make a cvterm, Sorry\n" and return unless ($params{name});
      my %cvtparams = (doc => $ldoc,
		       name => $params{name},
		      );
      delete $params{name};

      if ($params{cv_id}) {
	$cvtparams{cv_id} = $params{cv_id};
	delete $params{cv_id}; 
      } elsif ($params{cv}) {
	$cvtparams{cv} = $params{cv};
	delete $params{cv}; 
      } else {
	print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
      }
      $params{cvterm_id} = create_ch_cvterm(%cvtparams);      
    }

    #create a pub element if necessary
    if ($params{pub}) {
      $params{pub_id} = create_ch_pub(doc => $ldoc,
				      uniquename => $params{pub},
				     );
      delete $params{pub};
    }

    #now set required rank to 0 if not provided
    $params{rank} = '0' unless $params{rank};

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fct_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    return $fct_el;
}

# CREATE cell_line_cvtermprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from  cell_line_cvtermprop type cv or XML::DOM cvterm element required 
#        Note: will default to making a cell_lineprop from 'cell_line_cvtermprop type' cv unless 
#              cvname is provided
# cvname - string optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_cell_line_cvtermprop {
    my %params = @_;
    $params{parentname} = 'cell_line_cvterm';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE cell_line_dbxref element
# params
# doc - XML::DOM::Document required
# cell_line_id - macro cell_line id or XML::DOM cell_line element optionaal to create freestanding cell_line_dbxref
# dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
# accession - string required unless dbxref_id provided
# db_id - macro db id or XML::DOM db element required unless dbxref_id provided
# db - string name of db
# version - string optional
# description - string optional
# is_current - string 't' or 'f' boolean default = 't' so don't pass unless
#              this shoud be changed
sub create_ch_cell_line_dbxref {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document
  
  my $fd_el = $ldoc->createElement('cell_line_dbxref');
  
  if ($params{cell_line_id}) {
    $fd_el->appendChild(_build_element($ldoc,'cell_line_id',$params{cell_line_id}));
    delete $params{cell_line_id};
  }
  
  my $ic;
  if (exists($params{is_current})) { #assign value to a var and then remove from params
    if ($params{is_current}) {
      $ic = $params{is_current};
    } else {
      $ic = 'false';
    }
    delete $params{is_current};
  }
  
  # create a dbxref element if necessary
  unless ($params{dbxref_id}) {
    print "WARNING - missing required parameters, NO GO.\n" and return unless 
      ($params{accession} and ($params{db_id} or $params{db}));
    if ($params{db}) {
      $params{db_id} = create_ch_db(doc => $ldoc,
				    name => $params{db},
				   );
      delete $params{db};
    }
    
    
    $params{dbxref_id} = create_ch_dbxref(%params);
  }

    $fd_el->appendChild(_build_element($ldoc,'dbxref_id',$params{dbxref_id})); #add dbxref element
  $fd_el->appendChild(_build_element($ldoc,'is_current',$ic)) if $ic;

    return $fd_el;
}

# CREATE cell_line_feature element
# params
# doc - XML::DOM::Document required
# organism_id - macro id for organism or XML::DOM organism element
# genus - string
# species - string
# NOTE: you can use the generic paramaters in the following cases:
#       1.  you are only building either a cell_line or feature element and not both
#       2.  or both cell_line and feature have the same organism
#       otherwise use the prefixed parameters
# WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
# cell_line_id - macro id for cell_line or XML::DOM cell_line element
# cell_uniquename - string cell_line uniquename
# cell_organism_id - macro id for organism or XML::DOM organism element to link to cell_line
# cell_genus
# cell_species
# feature_id - macro id for feature or XML::DOM feature element
# feat_uniquename
# feat_organism_id - macro id for organism or XML::DOM organism element to link to feature
# feat_genus
# feat_species
# feat_type_id
# feat_type
# pub_id - macro pub id or XML::DOM pub element - required unless puname provided
# pub - string uniquename for pub (note will have lookup so can't create new pub here)
sub create_ch_cell_line_feature {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  #have to think about the checks a bit

  my $lf_el = $ldoc->createElement('cell_line_feature');

  # deal with feature bits
  if ($params{feat_uniquename}) {
    print "ERROR -- you don't have required parameters to make a feature, NO GO!\n" and return
      unless (($params{organism_id} or $params{feat_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{feat_genus} and $params{feat_species}))
	      and ($params{feat_type_id} or $params{feat_type}));

    unless ($params{feat_organism_id} or ($params{feat_genus} and $params{feat_species})) {
      $params{feat_organism_id} = $params{organism_id} if $params{organism_id};
      $params{feat_genus} = $params{genus} if $params{genus};
      $params{feat_species} = $params{species} if $params{species};
    }

    # gather all the feature parameters
    my %fparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /feat_(.+)/) {
	$fparams{$1} = $params{$p};
      }
    }

    $params{feature_id} = create_ch_feature(%fparams);
  }

  # likewise deal with the cell_line bits
  if ($params{cell_uniquename}) {
    print "ERROR -- you don't have required parameters to make a cell_line, NO GO!\n" and return
      unless (($params{organism_id} or $params{cell_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{cell_genus} and $params{cell_species})));

    unless ($params{cell_organism_id} or ($params{cell_genus} and $params{cell_species})) {
      $params{cell_organism_id} = $params{organism_id} if $params{organism_id};
      $params{cell_genus} = $params{genus} if $params{genus};
      $params{cell_species} = $params{species} if $params{species};
    }

    # gather all the cell_line parameters
    my %lparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /cell_(.+)/) {
	$lparams{$1} = $params{$p};
      }
    }

    $params{cell_line_id} = create_ch_cell_line(%lparams);
  }
  # finally add the pub info
  unless ($params{pub_id}) {
    print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
    $params{pub_id} = create_ch_pub(doc => $ldoc,
				    uniquename => $params{pub},
				   );
  }
  # and then add the feature, cell_line or both to the cell_line_feature element
  $lf_el->appendChild(_build_element($ldoc,'cell_line_id',$params{cell_line_id})) if $params{cell_line_id};
  $lf_el->appendChild(_build_element($ldoc,'feature_id',$params{feature_id})) if $params{feature_id};
  $lf_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id})) if $params{pub_id};
  return $lf_el;
}

# CREATE cell_line_library element
# params
# doc - XML::DOM::Document required
# organism_id - macro id for organism or XML::DOM organism element
# genus - string
# species - string
# NOTE: you can use the generic paramaters in the following cases:
#       1.  you are only building either a cell_line or library element and not both
#       2.  or both cell_line and library have the same organism
#       otherwise use the prefixed parameters
# WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
# cell_line_id - macro id for cell_line or XML::DOM cell_line element
# cell_uniquename - string cell_line uniquename
# cell_organism_id - macro id for organism or XML::DOM organism element to link to cell_line
# cell_genus
# cell_species
# library_id - macro id for library or XML::DOM library element
# lib_uniquename
# lib_organism_id - macro id for organism or XML::DOM organism element to link to library
# lib_genus
# lib_species
# lib_type_id
# lib_type
# pub_id - macro pub id or XML::DOM pub element - required unless puname provided
# pub - string uniquename for pub (note will have lookup so can't create new pub here)
sub create_ch_cell_line_library {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  #have to think about the checks a bit

  my $lf_el = $ldoc->createElement('cell_line_library');

  # deal with library bits
  if ($params{lib_uniquename}) {
    print "ERROR -- you don't have required parameters to make a library, NO GO!\n" and return
      unless (($params{organism_id} or $params{lib_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{lib_genus} and $params{lib_species}))
	      and ($params{lib_type_id} or $params{lib_type}));

    unless ($params{lib_organism_id} or ($params{lib_genus} and $params{lib_species})) {
      $params{lib_organism_id} = $params{organism_id} if $params{organism_id};
      $params{lib_genus} = $params{genus} if $params{genus};
      $params{lib_species} = $params{species} if $params{species};
    }

    # gather all the library parameters
    my %fparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /lib_(.+)/) {
	$fparams{$1} = $params{$p};
      }
    }

    $params{library_id} = create_ch_library(%fparams);
  }

  # likewise deal with the cell_line bits
  if ($params{cell_uniquename}) {
    print "ERROR -- you don't have required parameters to make a cell_line, NO GO!\n" and return
      unless (($params{organism_id} or $params{cell_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{cell_genus} and $params{cell_species})));

    unless ($params{cell_organism_id} or ($params{cell_genus} and $params{cell_species})) {
      $params{cell_organism_id} = $params{organism_id} if $params{organism_id};
      $params{cell_genus} = $params{genus} if $params{genus};
      $params{cell_species} = $params{species} if $params{species};
    }

    # gather all the cell_line parameters
    my %lparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /cell_(.+)/) {
	$lparams{$1} = $params{$p};
      }
    }

    $params{cell_line_id} = create_ch_cell_line(%lparams);
  }
  # finally add the pub info

  unless ($params{pub_id}) {
    print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
    $params{pub_id} = create_ch_pub(doc => $ldoc,
				    uniquename => $params{pub},
				   );
  }
  # and then add the library, cell_line or both to the cell_line_library element
  $lf_el->appendChild(_build_element($ldoc,'cell_line_id',$params{cell_line_id})) if $params{cell_line_id};
  $lf_el->appendChild(_build_element($ldoc,'library_id',$params{library_id})) if $params{library_id};
  $lf_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id})) if $params{pub_id};

  return $lf_el;
}

# CREATE cell_line_libraryprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for cell_line_libraryprop type or XML::DOM cvterm element required
# type -  string from cell_line_libraryprop type cv 
#        Note: will default to making a cell_line_library from 'cell_line_libraryprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'cell_line_libraryprop type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_cell_line_libraryprop {
    my %params = @_;
    $params{parentname} = 'cell_line_library';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}


# CREATE cell_line_pub element
# Note that this is just calling create_ch_pub  
#      and adding returned pub_id element as a child of cell_line_pub
#      or just appending the pub element to cell_line_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
#        creating a new pub (i.e. not null value but not part of unique key
# type - string for type from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_cell_line_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('cell_line_pub');
    if ($params{cell_line_id}) {
	$fp_el->appendChild(_build_element($ldoc,'cell_line_id',$params{cell_line_id}));
	delete $params{cell_line_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE cell_line_relationship element
# NOTE: this will create either a subject_id or object_id cell_line_relationship element
# but you have to attach this to the related cell_line elsewhere
# params
# doc - XML::DOM::Document required
# object _id - macro id for object cell_line or XML::DOM cell_line element
# subject_id - macro id for subject cell_line or XML::DOM cell_line element
# NOTE you can pass one or both of the above parameters with the following rules:
# if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other cell_line
# if both are passed then is_object, is_subject and any parameters to create a cell_line are ignored
# is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
# is_subject - boolean 't'          this flag indicates if the cell_line info provided should be 
#                                   added in as subject or object cell_line
# type_id - macro id for relationship type or XML::DOM cvterm element (Note: currently all is_relationship = '0'
# type - string for relationship type note: with this param  type will be assigned to relationship_type cv
# rank - integer optional with default = 0
# cell_line_id - macro_id for a cell_line or XML::DOM cell_line element required unless minimal cell_line bits provided
# uniquename - string required unless cell_line provided
# organism_id - macro id for organism or XML::DOM organism element required unless cell_line or (genus & species) provided
# genus - string required unless cell_line or organism provided
# species - string required unless cell_line or organism provided
sub create_ch_cell_line_relationship {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    ## cell_line_relationship element (will be returned)
    $params{relation} = $ldoc->createElement('cell_line_relationship');

    # need to deal with different param names for relationship type info
    if ($params{type}) {
      $params{rtype} = $params{type};
      delete $params{type};
    }
    if ($params{type_id}) {
      $params{rtype_id} = $params{type_id};
      delete $params{type_id};
    }
        
    # deal with creating the cell_line bits if present
    unless ($params{object_id} and $params{subject_id}) {
      unless ($params{cell_line_id}) {
	unless ($params{organism_id}) {
	_add_orgid_param(\%params);
	}
      	# before creating cell_line figure out which parameters
	my %fparams = (doc => $ldoc,
		       uniquename => $params{uniquename},
		       organism_id => $params{organism_id},
		      );
	delete $params{uniquename};
	delete $params{organism_id};
	$params{cell_line_id} = create_ch_cell_line(%fparams);
      } # now we have a cell_line element
      $params{thingy} = $params{cell_line_id};
      delete $params{cell_line_id};
    }
    $params{rtypecv} = 'cell_line_relationship' if ($params{rtype} and ! $params{rtypecv});   
    
    return create_ch_relationship(%params);
  }
*create_ch_clr = \&create_ch_cell_line_relationship;
*create_ch_c_l_r = \&create_ch_cell_line_relationship;

# CREATE cell_line_synonym element
# params
# doc - XML::DOM::Document required
# synonym_id - XML::DOM synonym element required unless name and type provided
# name - string required unless synonym_id element provided
# type_id - macro id for synonym type or XML::DOM cvterm element
#                 required unless a synonym element provided
# type - string = name from the 'synonym type' cv
# pub_id macro id for a pub or a XML::DOM pub element required
# pub - a pub uniquename (i.e. FBrf)
# synonym_sgml - string optional but if not provided then synonym_sgml = name
#               - do not provide if a synonym element is provided
# is_current - optional string = 'f' or 't' default is 't' so don't provide this param
#               unless you know you want to change the value
# is_internal - optional string = 't' default is 'f' so don't provide this param
#               unless you know you want to change the value
sub create_ch_cell_line_synonym {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $ls_el = $ldoc->createElement('cell_line_synonym');

    # this is exactly the same thing as create_ch_feature_synonym with different name
    # so call that method then change parentage?
    my $el = create_ch_feature_synonym(%params);
    my @children = $el->getChildNodes;
    $ls_el->appendChild($_) for @children;

    return $ls_el;
}

# CREATE contact or contact_id element
# params
# doc - XML::DOM::Document optional - required
# name - string required
# description - string optional
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value
# with_id - boolean optional if 1 then contact_id element is returned
sub create_ch_contact {
  my %params = @_;
  $params{elname} = 'contact';
  $params{required} = ['name'];
  my $eel = _create_simple_element(%params);
  return $eel;
}   

# CREATE cv or cv_id element
# params
# doc - XML::DOM::Document required
# name - string required
# definition - string optional
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value
# with_id - boolean optional if 1 then cv_id element is returned
sub create_ch_cv {
  my %params = @_;
  $params{elname} = 'cv';
  $params{required} = ['name'];
  my $eel = _create_simple_element(%params);
  return $eel;
}

# CREATE cvterm element
# params
# doc - XML::DOM::Document required
# name - string required
# cv_id - string = macro_id or XML::DOM cv element required
# cv - name of a cv
# definition - string optional
# dbxref_id - macro_id string or XML::DOM dbxref element
# is_obsolete - boolean optional default = false
# is_relationshiptype - boolean optional
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value 
# no_lookup - boolean option if 1 then default op="lookup" attribute will not be added to element
# note that we don't have a with_id parameter because either it will be freestanding 
# term or will have another type of id (e.g. type_id)
sub create_ch_cvterm {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    # check for required parameters
    print "WARNING -- missing parameters required for at least one of the two unique keys on cvterm\n"
      and return
	unless ($params{dbxref_id} or ($params{name} and ($params{cv_id} or $params{cv})));

    ## cvterm element (will be returned)
    my $cvt_el = $ldoc->createElement('cvterm');
    $cvt_el->setAttribute('id',$params{macro_id}) if $params{macro_id};  

    # add an op="lookup" attribute unless no_lookup is specified
    unless ($params{no_lookup}) {
	$cvt_el->setAttribute('op','lookup');
    }

    # check for cv parameter and if present convert into cv_id element
    if ($params{cv}) {
      $params{cv_id} = create_ch_cv(doc => $ldoc,
				    name => $params{cv},);
      delete $params{cv};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc' || $e eq 'macro_id' || $e eq 'no_lookup');
	$cvt_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
  
    return $cvt_el;
}

sub create_ch_cvtermprop {
    my %params = @_;
    $params{parentname} = 'cvterm';
    unless ($params{type_id}) {
	$params{cvname} = 'cvterm_property_type' unless $params{cvname};
    }
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}



# CREATE cvterm_relationship element
# NOTE: this will create either a subject_id or object_id cvterm_relationship but 
# you have to attach this to the related cvterm elsewhere.
# params
# doc - XML::DOM::Document required
# is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
# is_subject - boolean 't'          this flag indicates if the cvterm info provided should be 
#                                   added in as subject or object cvterm
# rtype_id - required macro id for relationship type cvterm or XML::DOM cvterm element unless rtype
# rtype - required string unless rtype_id provided
# note: if rtype is used cvterm will be assigned to 'relationship type' cv
# cvterm_id - macro id for cvterm or XML::DOM cvterm element required unless cvterm bits provided
# name - name of the cvterm - required unless cvterm element provided
# cv_id - macro id for cv or XML::DOM cv element required unless cv provided if you have a name
# cv - string name of a cv required unless cv_id if name provided
# dbxref_id - macro id for cvterm dbxref or XML::DOM dbxref element optional
# macro_id - string optional will add id attribute to cvterm element
sub create_ch_cvterm_relationship {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    ## cvterm_relationship element (will be returned)
  $params{relation} = $ldoc->createElement('cvterm_relationship');
    
    # deal with creating the cvterm bits if present
    unless ($params{object_id} and $params{subject_id}) {
      unless ($params{cvterm_id}) {
	print "WARNING -- missing required parameters for cvterm creation\n" and return
	  unless $params{name} and ($params{cv} or $params{cv_id});
	
	my %cvt_params = (doc => $ldoc,
			  name => $params{name},
			 );
	delete $params{name};
	
	if ($params{cv_id}) {
	  $cvt_params{cv_id} = $params{cv_id};
	  delete $params{cv_id};
	} else {
	  $cvt_params{cv} = $params{cv};
	  delete $params{cv};
	}
	if ($params{dbxref_id}) {
	  $cvt_params{dbxref_id} = $params{dbxref_id} ;
	  delete $params{dbxref_id};
	}
	if ($params{macro_id}){
	  $cvt_params{macro_id} = $params{macro_id};
	  delete $params{macro_id};
	}
	$params{cvterm_id} = create_ch_cvterm(%cvt_params);
      } # now we have a cvterm element to associate as subject or object
      $params{thingy} = $params{cvterm_id};
      delete $params{cvterm_id};
    }
    # not sure if there is a default cv for cvterm relationship types need to pass the rtypecv value?
    # for now it will default to 'relationship type'

    return create_ch_relationship(%params);
}
*create_ch_cr = \&create_ch_cvterm_relationship;

# CREATE db or db_id element
# params
# doc - XML::DOM::Document required
# name - string required
# contact_id - string = a macro id or XML::DOM contact element optional
# contact - string = contact name optional NOTE: if you provide both contact_id and contact the contact 
#                                                value will be used
# description - string optional
# urlprefix - string optional
# url - string optional
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value
# with_id - boolean optional if 1 then db_id element is returned
sub create_ch_db {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $dbid_el = $ldoc->createElement('db_id') if $params{with_id};
    my $db_el = $ldoc->createElement('db');
    $db_el->setAttribute('id',$params{macro_id}) if $params{macro_id};

    # check to see if contact param is used and if it is create a contact_id element
    # and remove contact from the param list
    if ($params{contact}) {
      $params{contact_id} = create_ch_contact(doc => $ldoc,
					      name => $params{contact},);
      delete $params{contact};
    }

    foreach my $e (keys %params) {
      next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id');
      $db_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }

    if ($dbid_el) {
      $dbid_el->appendChild($db_el);
      return $dbid_el;
    }
    return $db_el;
}


# CREATE dbxref or dbxref_id element
# params
# doc - XML::DOM::Document required
# accession - string required
# db_id - string = macro id for db or XML::DOM db element required
# db - string dbname required if db_id not provided
# version - string optional
# description - string optional
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value
# with_id - boolean optional if 1 then dbxref_id element is returned
# no_lookup - boolean option if 1 then default op="lookup" attribute will not be added to element
sub create_ch_dbxref {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    ## dbxref_id element (will be returned)
    my $dxid_el = $ldoc->createElement('dbxref_id') if $params{with_id};
    my $dx_el = $ldoc->createElement('dbxref');
    $dx_el->setAttribute('id',$params{macro_id}) if $params{macro_id};

    # add an op="lookup" attribute unless no_lookup is specified
    unless ($params{no_lookup}) {
	$dx_el->setAttribute('op','lookup');
    }

    # check for db param and if so make db_id element and delete db param
    if ($params{db}) {
      $params{db_id} = create_ch_db(doc => $ldoc,
				    name => $params{db},);
      delete $params{db};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id' || $e eq 'no_lookup');
	$dx_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }

    if ($dxid_el) {
	$dxid_el->appendChild($dx_el);
	return $dxid_el;
    }
    return $dx_el;
}

# CREATE dbxrefprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id or XML::DOM cvterm element required
#        Note: will default to making a dbxrefprop from 'property type' cv unless cvname is provided
# type - string from property type cv
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_dbxrefprop {
    my %params = @_;
    $params{parentname} = 'dbxref';
    unless ($params{type_id}) {
        $params{cvname} = 'property type' unless $params{cvname};
    }
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}


# CREATE environment or environment_id element
# params
# doc - XML::DOM::Document optional - required
# uniquename - string required
# description - string optional
# with_id - boolean optional if 1 then expression_id element is returned
# macro_id - string to specify as id attribute of this element for later use
sub create_ch_environment {
    my %params = @_;
    $params{elname} = 'environment';
    $params{required} = ['uniquename'];
    my $eel = _create_simple_element(%params);
    return $eel;
}

# CREATE environment_cvterm element
# params
# doc - XML::DOM::Document optional - required
# environment_id - macro id for environment of XML::DOM environment element
# uniquename - environment uniquename
# NOTE: you need to pass environment bits if attaching to existing cvterm element or 
#       creating a freestanding environment_cvterm
# cvterm_id -  macro id for cvterm of XML::DOM cvterm element
# name - cvterm name
# cv_id - macro id for a CV or XML::DOM cv element
# cv - name of a cv
# is_obsolete - optional param for cvterm
# NOTE: you need to pass cvterm bits if attaching to existing environment element or 
#       creating a freestanding environment_cvterm
sub create_ch_environment_cvterm {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  my $ect_el = $ldoc->createElement('environment_cvterm');

  # create an environment term if necessary
  if ($params{uniquename}) {
    $params{environment_id} = create_ch_environment(doc => $ldoc,
						    uniquename => $params{uniquename},
						   );
    delete $params{uniquename};
  }
    
  # create a cvterm element if necessary
  if ($params{name}) {
    print "ERROR: You don't have all the parameters required to make a cvterm, NO GO!\n" and return
      unless ($params{cv_id} or $params{cv});
    my %cvtparams = (doc => $ldoc,
		     name => $params{name},
		    );
    delete $params{name};
    
    if ($params{cv_id}) {
      $cvtparams{cv_id} = $params{cv_id};
      delete $params{cv_id}; 
    } elsif ($params{cv}) {
      $cvtparams{cv} = $params{cv};
      delete $params{cv}; 
    } else {
      print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
    }

    if ($params{is_obsolete}) {
      $cvtparams{is_obsolete} = $params{is_obsolete};
      delete $params{is_obsolete};
    }
    $params{cvterm_id} = create_ch_cvterm(%cvtparams);      
  }

  # now see which elements to attach to environment_cvterm
  $ect_el->appendChild(_build_element($ldoc,'environment_id',$params{environment_id})) if $params{environment_id};
  $ect_el->appendChild(_build_element($ldoc,'cvterm_id',$params{cvterm_id})) if $params{cvterm_id};
  
  return $ect_el;
}

# CREATE expression or expression_id element
# params
# doc - XML::DOM::Document optional - required
# uniquename - string required
# md5checksum
# description - string optional
# with_id - boolean optional if 1 then expression_id element is returned
# macro_id - string to specify as id attribute of this element for later use
sub create_ch_expression {
    my %params = @_;
    $params{elname} = 'expression';
    $params{required} = ['uniquename'];
    my $eel = _create_simple_element(%params);
    return $eel;
}

# CREATE expression_cvterm element
# params
# doc - XML::DOM::Document required
# cvterm_id - XML::DOM cvterm element unless other cvterm bits are provided
# name - string required unless cvterm_id provided
# cv_id -  macro id for a cv or XML::DOM cv element required unless cvterm_id or cv provided
# cv - string = cvname required unless cvterm_id or cv_id provided
# cvterm_type_id - macro id for expression slots cvterm or XML::DOM cvterm element - required unless cvterm_type
# cvterm_type - string from the expression slots cv
# rank - int (default = 0 so don't pass unless you know it's different)  
# think about adding is_not column with default = false
sub create_ch_expression_cvterm {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  my $ect_el = $ldoc->createElement('expression_cvterm');

  # create a cvterm element if necessary
  unless ($params{cvterm_id}) {
    print "WARNING -- you haven't provided params to make a cvterm, Sorry\n" and return unless ($params{name});
    my %cvtparams = (doc => $ldoc,
		     name => $params{name},
		    );
    delete $params{name};
    
    if ($params{cv_id}) {
      $cvtparams{cv_id} = $params{cv_id};
      delete $params{cv_id}; 
    } elsif ($params{cv}) {
      $cvtparams{cv} = $params{cv};
      delete $params{cv}; 
    } else {
      print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
    }
    $params{cvterm_id} = create_ch_cvterm(%cvtparams);      
  }
  
  # deal with the type_id info 
  unless ($params{cvterm_type_id}) {
      print "WARNING -- you are trying to create a expression type without providing the info - Sorry, NO GO\n"
	and return unless $params{cvterm_type};

      $params{cvterm_type_id} = create_ch_cvterm(doc => $ldoc,
						 name => $params{cvterm_type},
						 cv => 'expression slots', # need to determine what this cv will be called
						);

      delete $params{cvterm_type};
  }
  
  #now set required rank to 0 if not provided
  $params{rank} = '0' unless $params{rank};

  foreach my $e (keys %params) {
    next if ($e eq 'doc');
    $ect_el->appendChild(_build_element($ldoc,$e,$params{$e}));
  }
  return $ect_el;
}

# CREATE expression_cvtermprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for cvterm from expression_cvtprop type cv or XML::DOM cvterm element required
# type -  string from  expression_cvtermprop type cv
#         Note: will default to making a expressionprop from 'expression_cvtermprop type' cv unless
#              cvname is provided
# cvname - string optional but see above for type and do not provide if passing type_id
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_expression_cvtermprop {
    my %params = @_;
    $params{parentname} = 'expression_cvterm';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE expression_pub element
# Note that this is just calling create_ch_pub setting with_id = 1
#      and adding returned pub_id element as a child of expression_pub
#      or just appending the pub element to expression_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
#            creating a new pub (i.e. not null value but not part of unique key
# type - string from pub type cv same requirement as type_id
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_expression_pub {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document
  
  my $ep_el = $ldoc->createElement('expression_pub');
  
    if ($params{expression_id}) {
	$ep_el->appendChild(_build_element($ldoc,'feature_id',$params{expression_id}));
	delete $params{expression_id};
    }

  unless ($params{pub_id}) {
    $params{pub_id} = create_ch_pub(%params); #will return a pub element
  }
  $ep_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));

  return $ep_el;
}

# CREATE expressionprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from realtionship property type cv or XML::DOM cvterm element required
#        Note: will default to making a expressionprop from 'property type' cv unless cvname is provided
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_expressionprop {
    my %params = @_;
    $params{parentname} = 'expression';
    unless ($params{type_id}) {
	$params{cvname} = 'expressionprop type' unless $params{cvname};
    }
    my $ep_el = create_ch_prop(%params);
    return $ep_el;
}

# CREATE feature element
# params
# doc - XML::DOM::Document required
# uniquename - string required
# type_id - cvterm macro id or XML::DOM cvterm element required unless type
# type - string for a cvterm Note: will default to using SO cv unless a cvname is provided
# cvname - string optional to specify a cv other than SO for the type_id
#          do not use if providing a cvterm element or macro id to type_id
# organism_id - organism macro id or XML::DOM organism element required if no genus and species
# genus - string required if no organism
# species - string required if no organism
# dbxref_id - dbxref macro id or XML::DOM dbxref element optional
# name - string optional
# residue - string optional
# seqlen - integer optional (if seqlen = 0 pass as string)
# md5checksum - string optional
# is_analysis - boolean 't' or 'f' default = 'f' optional
# is_obsolete - boolean 't' or 'f' default = 'f' optional
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value
# with_id - boolean optional if 1 then feature_id element is returned
# no_lookup - boolean option if 1 then default op="lookup" attribute will not be added to element
sub create_ch_feature {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fid_el = $ldoc->createElement('feature_id') if $params{with_id};

    ## feature element (will be returned)
    my $f_el = $ldoc->createElement('feature');
    $f_el->setAttribute('id',$params{macro_id}) if $params{macro_id};    

    # add an op="lookup" attribute unless no_lookup is specified
    unless ($params{no_lookup}) {
	$f_el->setAttribute('op','lookup');
    }
	
    #create organism_id element if genus and species are provided
    unless ($params{organism_id}) {
	_add_orgid_param(\%params);
    }

    # figure out which cv to use in type_id element that we make below
    my $cv = 'SO';
    if ($params{cvname}) {
	$cv = $params{cvname};
	delete $params{cvname};
    }

    # now deal make a cvterm element for type_id if string is provided
    if ($params{type}) {
      $params{type_id} = create_ch_cvterm(doc => $ldoc,
				   name => $params{type},
				   cv => $cv,
				  );
      delete $params{type};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id' || $e eq 'no_lookup');
	$f_el->appendChild(_build_element($ldoc, $e,$params{$e}));
    }

    if ($fid_el) {
	$fid_el->appendChild($f_el);
	return $fid_el;
    }
    return $f_el;
}


# CREATE feature_cvterm element
# params
# doc - XML::DOM::Document required
# cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
# name - string required unless cvterm_id provided Note: a cvterm has a lookup by default cannot make a new cvterm 
#                                                        with this method
# cv_id - macro id for cv or XML::DOM cv element required if name and not cv
# cv - string for name of cv required if name and not cv_id
# pub_id - macro id for pub or XML::DOM pub element required unless pub
# pub - string = pub uniquename Note: as pub has lookup option by default can't make a new pub using this param
# is_not - optional boolean 't' or 'f' with default = 'f' so don't pass unless you know you want to change
sub create_ch_feature_cvterm {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fct_el = $ldoc->createElement('feature_cvterm');

    #create a cvterm element if necessary
    unless ($params{cvterm_id}) {
      print "WARNING -- you haven't provided params to make a cvterm, Sorry\n" and return unless ($params{name});
      my %cvtparams = (doc => $ldoc,
		       name => $params{name},
		      );
      delete $params{name};

      if ($params{cv_id}) {
	$cvtparams{cv_id} = $params{cv_id};
	delete $params{cv_id}; 
      } elsif ($params{cv}) {
	$cvtparams{cv} = $params{cv};
	delete $params{cv}; 
      } else {
	print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
      }
      $params{cvterm_id} = create_ch_cvterm(%cvtparams);      
    }

    #create a pub element if necessary
    if ($params{pub}) {
      $params{pub_id} = create_ch_pub(doc => $ldoc,
				      uniquename => $params{pub},
				     );
      delete $params{pub};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fct_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    return $fct_el;
}

# CREATE feature_cvtermprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from  feature_cvtermprop type cv or XML::DOM cvterm element required 
#        Note: will default to making a featureprop from 'feature_cvtermprop type' cv unless 
#              cvname is provided
# cvname - string optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_feature_cvtermprop {
    my %params = @_;
    $params{parentname} = 'feature_cvterm';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}


# CREATE feature_dbxref element
# params
# doc - XML::DOM::Document required
# feature_id - macro feature id or XML::DOM feature element optionaal to create freestanding feature_dbxref
# dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
# accession - string required unless dbxref_id provided
# db_id - macro db id or XML::DOM db element required unless dbxref_id provided
# db - string name of db
# version - string optional
# description - string optional
# is_current - string 't' or 'f' boolean default = 't' so don't pass unless
#              this shoud be changed
sub create_ch_feature_dbxref {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fd_el = $ldoc->createElement('feature_dbxref');

    if ($params{feature_id}) {
	$fd_el->appendChild(_build_element($ldoc,'feature_id',$params{feature_id}));
	delete $params{feature_id};
    }

    my $ic;
    if (exists($params{is_current})) { #assign value to a var and then remove from params
      if ($params{is_current}) {
	$ic = $params{is_current};
      } else {
	$ic = 'false';
      }
      delete $params{is_current};
    }
    
    # create a dbxref element if necessary
    unless ($params{dbxref_id}) {
      print "WARNING - missing required parameters, NO GO.\n" and return unless 
	($params{accession} and ($params{db_id} or $params{db}));
      if ($params{db}) {
	$params{db_id} = create_ch_db(doc => $ldoc,
				      name => $params{db},
				     );
	delete $params{db};
      }
      
	
      $params{dbxref_id} = create_ch_dbxref(%params);
    }

    $fd_el->appendChild(_build_element($ldoc,'dbxref_id',$params{dbxref_id})); #add dbxref element
    $fd_el->appendChild(_build_element($ldoc,'is_current',$ic)) if $ic;

    return $fd_el;
}

# CREATE feature_expression element
# params
# doc - XML::DOM::Document required
# feature_id - OPTIONAL macro feature id or XML::DOM feature element to create freestanding feature_expression 
# expression_id - macro expression id or XML::DOM expression element - required unless uniquename provided
# uniquename - string required unless expression_id
# pub_id - macro pub id or XML::DOM pub element - required unless puname provided
# pub - string uniquename for pub (note will have lookup so can't create new pub here)
sub create_ch_feature_expression {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fd_el = $ldoc->createElement('feature_expression');

    # create a expression element if necessary
    unless ($params{expression_id}) {
      print "WARNING - missing required expression info, NO GO.\n" and return unless $params{uniquename};
      $params{expression_id} = create_ch_expression(doc => $ldoc,
						   uniquename => $params{uniquename},);
      delete $params{uniquename};
    }

    #create a pub element if necessary
    unless ($params{pub_id}) {
	print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
	$params{pub_id} = create_ch_pub(doc => $ldoc,
					uniquename => $params{pub},
				       );
	delete $params{pub};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fd_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }

    return $fd_el;
}

# CREATE feature_expressionprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for feature_expressionprop type or XML::DOM cvterm element required
# type -  string from  expression_cvtermprop type cv 
#        Note: will default to making a featureprop from 'expression_cvtermprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'feature_expressionprop type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_feature_expressionprop {
    my %params = @_;
    $params{parentname} = 'feature_expression';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE feature_genotype element
# params
# doc - XML::DOM::Document optional - required
# genotype_id - macro_id for a genotype or XML::DOM genotype element
# uniquename - string = genotype uniquename
# chromosome_id macro id for chromosome feature or XML::DOM feature element
#      NOTE: that this is part of the key but can be null so not required
# cgroup - int (default = 0 so don't pass unless you know its different)
# rank - int (default = 0 so don't pass unless you know its different)
# cvterm_id - macro id for a cvterm or XML::DOM cvterm element required
sub create_ch_feature_genotype {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    ## feature_genotype element (will be returned)
    my $fl_el = $ldoc->createElement('feature_genotype');

    #first warn about potential oddness
    print STDERR "WARNING - no genotype_id object\n" and return unless ($params{genotype_id} or $params{uniquename});
    print STDERR "WARNING - no chromosome_id feature object\n" unless $params{chromosome_id};
    print STDERR "WARNING - no cvterm_id object\n" and return unless $params{cvterm_id};

    #now set required cgroup and ranks to 0 if not provided
    $params{cgroup} = '0' unless $params{cgroup};
    $params{rank} = '0' unless $params{rank};

    # make a genotype element if uniquename is specified
    if ($params{uniquename}) {
      $params{genotype_id} = create_ch_genotype(doc => $ldoc,
						uniquename => $params{uniquename},
					       );
      delete $params{uniquename};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fl_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    return $fl_el
}

# CREATE feature_grpmember element
# params 
# doc - XML::DOM::Document optional - required
#
# Here are parameters to make a grpmember element must either pass a grpmember_id or the other necessary bits
# grpmember_id - macro grpmember id or XML::DOM grpmember element
#
# Here are parameters to make an feature element must either pass feature_id or required bits
# feature_id - macro feature id or XML::DOM feature element

sub create_ch_feature_grpmember {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $af_el = $ldoc->createElement('feature_grpmember');
     
    # or if a macro id or existing element have been provided
    if ($params{grpmember_id}) {
      $af_el->appendChild(_build_element($ldoc,'grpmember_id',$params{grpmember_id}));
      delete $params{grpmember_id};
    }

   # and here we are dealing with feature info if provided
    if ($params{feature_id}) {
      $af_el->appendChild(_build_element($ldoc,'feature_id',$params{feature_id}));
      delete $params{feature_id};
    }

    return $af_el;
}
# CREATE feature_grpmember_pub element
# Note that this is just calling create_ch_pub setting with_id = 1 
#      and adding returned pub_id element as a child of feature_grpmember_pub
#      or just appending the pub element to feature_grpmember_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - macro_id for pub or XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
#        creating a new pub (i.e. not null value but not part of unique key
# type - string for type from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_feature_grpmember_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('feature_grpmember_pub');

    if ($params{feature_grpmember_id}) {
	$fp_el->appendChild(_build_element($ldoc,'feature_grpmember_id',$params{feature_grpmember_id}));
	delete $params{feature_grpmember_id};
    }
    
    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE feature_humanhealth_dbxref element
# params
# doc - XML::DOM::Document required
# feature_id - OPTIONAL macro feature id or XML::DOM feature element to create freestanding feature_humanhealth_dbxref 
# humanhealth_dbxref_id - macro humanhealth_dbxref id or XML::DOM humanhealth_dbxref element - required
# pub_id - macro pub id or XML::DOM pub element - required unless puname provided
# pub - string uniquename for pub (note will have lookup so can't create new pub here)
sub create_ch_feature_humanhealth_dbxref {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fd_el = $ldoc->createElement('feature_humanhealth_dbxref');

    #first warn about potential oddness
    print STDERR "WARNING - no pub_id object\n" and return unless ($params{pub_id} or $params{pub});
    print STDERR "WARNING - no feature_id object\n" unless $params{feature_id};
    print STDERR "WARNING - no humanhealth_dbxref_id object\n" and return unless $params{humanhealth_dbxref_id};

    #create a pub element if necessary
    unless ($params{pub_id}) {
	print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
	$params{pub_id} = create_ch_pub(doc => $ldoc,
					uniquename => $params{pub},
				       );
	delete $params{pub};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fd_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }

    return $fd_el;
}

# CREATE feature_interaction element
# params
# doc - XML::DOM::Document required
# feature_id - OPTIONAL macro feature id or XML::DOM feature element to create freestanding feature_expression 
# interaction_id - macro interaction_id or XML::DOM interaction element - required unless uniquename and type info provided
# uniquename - string required unless interaction_id
# type_id - macro type_id or cvterm  element required unless interaction_id or type
# type - string required unless interaction_id or type_id
# cvname - string optional
# role_id - macro role_id or cvterm element unless role
# role - string term from 'PSI-MI' cv
# rank - int optional default = 0
sub create_ch_feature_interaction {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fi_el = $ldoc->createElement('feature_interaction');

    # create a interaction element if necessary
    unless ($params{interaction_id}) {
      print "WARNING - missing required expression info, NO GO.\n" and return 
	unless ($params{uniquename} and ($params{type_id} or $params{type}));

      # deal with type
      # now deal make a cvterm element for type_id if string is provided
      if ($params{type}) {
	my $cv = 'PSI-MI';
	if ($params{cvname}) {
	  $cv = $params{cvname};
	  delete $params{cvname};
	}
	$params{type_id} = create_ch_cvterm(doc => $ldoc,
					    name => $params{type},
					    cv => $cv,
					   );
	delete $params{type};
      }
      $params{interaction_id} = create_ch_interaction(doc => $ldoc,
						      uniquename => $params{uniquename},
						      type_id => $params{type_id},
						     );
      delete $params{uniquename};
      delete $params{type_id};
    }

    #create a role_id element if necessary
    unless ($params{role_id}) {
      print "WARNING - missing required role info, NO GO.\n" and return unless $params{role};
      $params{role_id} = create_ch_cvterm(doc => $ldoc,
					  name => $params{role},
					  cv => 'PSI-MI',
				       );
	delete $params{role};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
#	print "$e\t$params{$e}\n";
	my $eel = _build_element($ldoc,$e,$params{$e});
	$fi_el->appendChild($eel);
    }

    return $fi_el;
}

# CREATE feature_interactionprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from interaction property type cv or XML::DOM cvterm element required
#        Note: will default to making a from 'interaction property type' cv unless cvname is provided
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_feature_interactionprop {
    my %params = @_;
    $params{parentname} = 'feature_interaction';
    unless ($params{type_id}) {
	$params{cvname} = 'feature_interaction property type' unless $params{cvname};
    }
    my $ep_el = create_ch_prop(%params);
    return $ep_el;
}



# CREATE feature_interaction_pub element
# Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of feature_interaction_pub
#      or just appending the pub element to feature_interaction_pub if that is passed
# params
# doc - XML::DOM::Document required
# feature_interaction_id - optional feature_interaction element or macro_id 
#          to create freestanding ele
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_feature_interaction_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('feature_interaction_pub');

    if ($params{feature_interaction_id}) {
	$fp_el->appendChild(_build_element($ldoc,'feature_interaction_id',$params{feature_interaction_id}));
	delete $params{feature_interaction_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE feature_pub element
# Note that this is just calling create_ch_pub setting with_id = 1 
#      and adding returned pub_id element as a child of feature_pub
#      or just appending the pub element to feature_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - macro_id for pub or XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
#        creating a new pub (i.e. not null value but not part of unique key
# type - string for type from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_feature_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('feature_pub');

    if ($params{feature_id}) {
	$fp_el->appendChild(_build_element($ldoc,'feature_id',$params{feature_id}));
	delete $params{feature_id};
    }
    
    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# just calling create_ch_prop with correct params to make desired prop
# note that there are cases here where the value is null
sub create_ch_feature_pubprop {
    my %params = @_;
    $params{parentname} = 'feature_pub';
    unless ($params{type_id}) {
	$params{cvname} = 'pubprop type' unless $params{cvname};
    }
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE feature_relationship element
# NOTE: this will create either a subject_id or object_id feature_relationship element
# but you have to attach this to the related feature elsewhere
# params
# doc - XML::DOM::Document required
# object _id - macro id for object feature or XML::DOM feature element
# subject_id - macro id for subject feature or XML::DOM feature element
# NOTE you can pass one or both of the above parameters with the following rules:
# if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other feature
# if both are passed then is_object, is_subject and any parameters to create a feature are ignored
# is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
# is_subject - boolean 't'          this flag indicates if the feature info provided should be 
#                                   added in as subject or object feature
# rtype_id - macro id for relationship type or XML::DOM cvterm element (Note: currently all is_relationship = '0'
# rtype - string for relationship type note: with this param  type will be assigned to relationship_type cv
# rank - integer optional with default = 0
# feature_id - macro_id for a feature or XML::DOM feature element required unless minimal feature bits provided
# uniquename - string required unless feature provided
# organism_id - macro id for organism or XML::DOM organism element required unless feature or (genus & species) provided
# genus - string required unless feature or organism provided
# species - string required unless feature or organism provided
# ftype_id - macro id for feature type or XML::DOM cvterm element required unless feature provided
# ftype - string for feature type (will be assigned to SO cv and can't be new as will be lookup)
sub create_ch_fr {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document
  
  ## feature_relationship element (will be returned)
  $params{relation} = $ldoc->createElement('feature_relationship');
  
  # deal with creating the strain bits if present and needed
  unless ($params{object_id} and $params{subject_id}) {
    # deal with creating the feature bits if present
    unless ($params{feature_id}) {
      unless ($params{organism_id}) {
	_add_orgid_param(\%params);
      }
      my %fparams = (doc => $ldoc,
		     uniquename => $params{uniquename},
		     organism_id => $params{organism_id},
		    );
      delete $params{uniquename};
      delete $params{organism_id};
      if ($params{ftype_id}) { 
	$fparams{type_id} = $params{ftype_id};
	delete $params{ftype_id}; 
      } elsif ($params{ftype}) {
	$fparams{type} = $params{ftype};
	delete $params{ftype};
      } else {
	print "WARNING -- you need to provide a feature type to make a feature!\n" and return;
	}
      $params{feature_id} = create_ch_feature(%fparams);
    } # now we have a feature element to associate as subject or object
    $params{thingy} = $params{feature_id};
    delete $params{feature_id};
  }
    return create_ch_relationship(%params);
}
*create_ch_feature_relationship = \&create_ch_fr;
*create_ch_f_r = \&create_ch_fr;

# CREATE feature_relationshipprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from feature_relationshipprop type cv or XML::DOM cvterm element required 
#        Note: will default to making a featureprop from 'fr property type' cv unless cvname is provided
# cvname - string optional 
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_frprop {
    my %params = @_;
    $params{parentname} = 'feature_relationship';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}
*create_ch_fr_prop = \&create_ch_frprop;
*create_ch_feature_relationshipprop = \&create_ch_frprop;
*create_ch_f_rprop = \&create_ch_frprop;

# CREATE feature_relationshipprop_pub element
# # Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of feature_relationshipprop_pub
#      or just appending the pub element to feature_relationshipprop_pub if that is passed
# params
# doc - XML::DOM::Document required
# feature_relationshipprop_id -optional feature_relationshipprop XML::DOM element or macro id
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_frprop_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('feature_relationshipprop_pub');

    if ($params{feature_relationshipprop_id}) {
	$fp_el->appendChild(_build_element($ldoc,'feature_relationshipprop_id',$params{feature_relationshipprop_id}));
	delete $params{feature_relationshipprop_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}
*create_ch_feature_relationshipprop_pub = \& create_ch_frprop_pub;
*create_ch_f_rprop_pub = \& create_ch_frprop_pub;

# CREATE feature_relationship_pub element
# # Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of feature_relationship_pub
#      or just appending the pub element to feature_relationship_pub if that is passed
# params
# doc - XML::DOM::Document required
# feature_relationship_id -optional feature_relationship XML::DOM element or macro id
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_fr_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('feature_relationship_pub');

    if ($params{feature_relationship_id}) {
	$fp_el->appendChild(_build_element($ldoc,'feature_relationship_id',$params{feature_relationship_id}));
	delete $params{feature_relationship_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}
*create_ch_feature_relationship_pub = \& create_ch_fr_pub;

# CREATE feature_synonym element
# params
# doc - XML::DOM::Document required
# synonym_id - macro id for synonym? or XML::DOM synonym element required unless name and type provided
# name - string required unless synonym_id element provided
# type_id - macro id for synonym type or XML::DOM cvterm element
# type - string = name from the 'synonym type' cv
# pub_id - macro id for pub or a XML::DOM pub element required
# pub - a pub uniquename (i.e. FBrf)
# synonym_sgml - string optional but if not provided then synonym_sgml = name
#               - do not provide if a synonym element is provided
# is_current - optional string = 'f' or 't' default is 't' so don't provide this param
#               unless you know you want to change the value
# is_internal - optional string = 't' default is 'f' so don't provide this param
#               unless you know you want to change the value
sub create_ch_feature_synonym {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    ## feature_synonym element (will be returned)
    my $fs_el = $ldoc->createElement('feature_synonym');

    #create a synonym element and undefine name, type and synonym_sgml bits
    unless($params{synonym_id}) {
      # gather params for synonym
      print "WARNING - you haven't provided info to create a synonym element\n" and return unless $params{name};
      my %syn_params = (doc => $ldoc,
			name => $params{name},
		       );
      delete $params{name};
      if ($params{synonym_sgml}) {
	$syn_params{synonym_sgml} = $params{synonym_sgml};
	delete $params{synonym_sgml};
      }
       
      # check for type or type_id
      if ($params{type}) {
	$syn_params{type} = $params{type};
	delete $params{type};
      } elsif ($params{type_id}) {
	$syn_params{type_id} = $params{type_id};
	delete $params{type_id};
      } else {
	print "WARNING - you haven't provided a synonym type\n" and return;
      }
      my $syn_el = create_ch_synonym(%syn_params);
      $params{synonym_id} = $syn_el;	
    }

    # check for pub 
    if ($params{pub}) {
      $params{pub_id} = create_ch_pub(doc => $ldoc,
				      uniquename => $params{pub},
				     );
      delete $params{pub};
    }

     print "WARNING - you need to provide pub info for the synonym link\n" and return unless $params{pub_id};

    foreach my $e (keys %params) {
      next if ($e eq 'doc');
      $fs_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    
    return $fs_el;
}


# CREATE featureloc element
# params
# none of these parameters are strictly required and some warning is done
# but if you misbehave you could set up some funky featurelocs?
# srcfeature_id macro id for a feature or XML::DOM feature element
# fmin - integer
# fmax - integer
# strand - 1, 1 or 0 (0 must be passed as string or else will be undef)
# phase - int
# residue_info - string
# locgroup - int (default = 0 so don't pass unless you know its different)
# rank - int (default = 0 so don't pass unless you know its different)
# is_fmin_partial - boolean 't' or 'f' default = 'f'
# is_fmax_partial - boolean 't' or 'f' default = 'f'
sub create_ch_featureloc {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    ## featureloc element (will be returned)
    my $fl_el = $ldoc->createElement('featureloc');

    #first warn about potential oddness
#    print STDERR "WARNING - no srcfeature object\n" unless $params{srcfeature_id};
#    print STDERR "WARNING - no min and max coordinate pair\n" unless ($params{fmin} && $params{fmax});
#    print STDERR "WARNING - no strand\n" unless $params{strand};

    #now set required loc_group and ranks to 0 if not provided
    $params{locgroup} = '0' unless $params{locgroup};
    $params{rank} = '0' unless $params{rank};

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fl_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    return $fl_el
}

# CREATE featureloc_pub element
# Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of featureloc_pub
#      or just appending the pub element to featureloc_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - macro id for pub or XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required unless pub_id
# type_id - macro id for pub or XML::DOM cvterm element optional unless 
#        creating a new pub (i.e. not null value but not part of unique key
# type - string from pub type cv
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_featureloc_pub {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document
  
  my $fp_el = $ldoc->createElement('featureloc_pub');
  # add optional featureloc_id element if param is passed
  if ($params{featureloc_id}) {
	$fp_el->appendChild(_build_element($ldoc,'featureloc_id',$params{featureloc_id}));
	delete $params{featureloc_id};
    }

  unless ($params{pub_id}) {
    $params{pub_id} = create_ch_pub(%params); #will return a pub element
  }
  $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
  return $fp_el;  
}

# CREATE featureprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id or XML::DOM cvterm element required 
#        Note: will default to making a featureprop from 'property type' cv unless cvname is provided
# type - string from property type cv
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_featureprop {
    my %params = @_;
    $params{parentname} = 'feature';
    unless ($params{type_id}) {
	$params{cvname} = 'property type' unless $params{cvname};
    }
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}


# CREATE featureprop_pub element
# Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of featureprop_pub
#      or just appending the pub element to featureprop_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_featureprop_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('featureprop_pub');

    if ($params{featureprop_id}) {
	$fp_el->appendChild(_build_element($ldoc,'featureprop_id',$params{featureprop_id}));
	delete $params{featureprop_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE genotype or genotype_id element
# params
# doc - XML::DOM::Document optional - required
# uniquename - string required
# description - string optional
# name - string optional
# macro_id
# with_id - boolean optional if 1 then genotype_id element is returned
sub create_ch_genotype {
  my %params = @_;
  $params{elname} = 'genotype';
  $params{required} = ['uniquename'];
  my $eel = _create_simple_element(%params);
  return $eel;
}

# CREATE grp element
# params
# doc - XML::DOM::Document required
# uniquename - string required
# type_id - cvterm macro id or XML::DOM cvterm element required unless type
# type - string for a cvterm Note: will default to using SO cv unless a cvname is provided
# cvname - string optional to specify a cv other than SO for the type_id
#          do not use if providing a cvterm element or macro id to type_id
# name - string optional
# is_analysis - boolean 't' or 'f' default = 'f' optional
# is_obsolete - boolean 't' or 'f' default = 'f' optional
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value
# with_id - boolean optional if 1 then grp_id element is returned
# no_lookup - boolean option if 1 then default op="lookup" attribute will not be added to element
sub create_ch_grp {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fid_el = $ldoc->createElement('grp_id') if $params{with_id};

    ## grp element (will be returned)
    my $f_el = $ldoc->createElement('grp');
    $f_el->setAttribute('id',$params{macro_id}) if $params{macro_id};    

    # add an op="lookup" attribute unless no_lookup is specified
    unless ($params{no_lookup}) {
	$f_el->setAttribute('op','lookup');
    }
	
    # figure out which cv to use in type_id element that we make below
    my $cv = 'SO';
    if ($params{cvname}) {
	$cv = $params{cvname};
	delete $params{cvname};
    }

    # now deal make a cvterm element for type_id if string is provided
    if ($params{type}) {
      $params{type_id} = create_ch_cvterm(doc => $ldoc,
				   name => $params{type},
				   cv => $cv,
				  );
      delete $params{type};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id' || $e eq 'no_lookup');
	$f_el->appendChild(_build_element($ldoc, $e,$params{$e}));
    }

    if ($fid_el) {
	$fid_el->appendChild($f_el);
	return $fid_el;
    }
    return $f_el;
}

# CREATE grp_cvterm element
# params
# doc - XML::DOM::Document required
# cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
# name - string required unless cvterm_id provided Note: a cvterm has a lookup by default cannot make a new cvterm 
#                                                        with this method
# cv_id - macro id for cv or XML::DOM cv element required if name and not cv
# cv - string for name of cv required if name and not cv_id
# is_not - optional boolean 't' or 'f' with default = 'f' so don't pass unless you know you want to change
sub create_ch_grp_cvterm {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fct_el = $ldoc->createElement('grp_cvterm');

    #create a cvterm element if necessary
    unless ($params{cvterm_id}) {
      print "WARNING -- you haven't provided params to make a cvterm, Sorry\n" and return unless ($params{name});
      my %cvtparams = (doc => $ldoc,
		       name => $params{name},
		      );
      delete $params{name};

      if ($params{cv_id}) {
	$cvtparams{cv_id} = $params{cv_id};
	delete $params{cv_id}; 
      } elsif ($params{cv}) {
	$cvtparams{cv} = $params{cv};
	delete $params{cv}; 
      } else {
	print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
      }
      $params{cvterm_id} = create_ch_cvterm(%cvtparams);      
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fct_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    return $fct_el;
}

# CREATE grp_dbxref element
# params
# doc - XML::DOM::Document required
# grp_id - macro grp id or XML::DOM grp element optionaal to create freestanding grp_dbxref
# dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
# accession - string required unless dbxref_id provided
# db_id - macro db id or XML::DOM db element required unless dbxref_id provided
# db - string name of db
# version - string optional
# description - string optional
# is_current - string 't' or 'f' boolean default = 't' so don't pass unless
#              this shoud be changed
sub create_ch_grp_dbxref {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fd_el = $ldoc->createElement('grp_dbxref');

    if ($params{grp_id}) {
	$fd_el->appendChild(_build_element($ldoc,'grp_id',$params{grp_id}));
	delete $params{grp_id};
    }

    my $ic;
    if (exists($params{is_current})) { #assign value to a var and then remove from params
      if ($params{is_current}) {
	$ic = $params{is_current};
      } else {
	$ic = 'false';
      }
      delete $params{is_current};
    }
    
    # create a dbxref element if necessary
    unless ($params{dbxref_id}) {
      print "WARNING - missing required parameters, NO GO.\n" and return unless 
	($params{accession} and ($params{db_id} or $params{db}));
      if ($params{db}) {
	$params{db_id} = create_ch_db(doc => $ldoc,
				      name => $params{db},
				     );
	delete $params{db};
      }
      
	
      $params{dbxref_id} = create_ch_dbxref(%params);
    }

    $fd_el->appendChild(_build_element($ldoc,'dbxref_id',$params{dbxref_id})); #add dbxref element
    $fd_el->appendChild(_build_element($ldoc,'is_current',$ic)) if $ic;

    return $fd_el;
}

# CREATE grpprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id or XML::DOM cvterm element required 
#        Note: will default to making a grpprop from 'property type' cv unless cvname is provided
# type - string from property type cv
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_grpprop {
    my %params = @_;
    $params{parentname} = 'grp';
    unless ($params{type_id}) {
	$params{cvname} = 'property type' unless $params{cvname};
    }
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}


# CREATE grpprop_pub element
# Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of grpprop_pub
#      or just appending the pub element to grpprop_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_grpprop_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('grpprop_pub');

    if ($params{grpprop_id}) {
	$fp_el->appendChild(_build_element($ldoc,'grpprop_id',$params{grpprop_id}));
	delete $params{grpprop_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE grp_pub element
# Note that this is just calling create_ch_pub setting with_id = 1 
#      and adding returned pub_id element as a child of grp_pub
#      or just appending the pub element to grp_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - macro_id for pub or XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
#        creating a new pub (i.e. not null value but not part of unique key
# type - string for type from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_grp_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('grp_pub');

    if ($params{grp_id}) {
	$fp_el->appendChild(_build_element($ldoc,'grp_id',$params{grp_id}));
	delete $params{grp_id};
    }
    
    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# just calling create_ch_prop with correct params to make desired prop
# note that there are cases here where the value is null
sub create_ch_grp_pubprop {
    my %params = @_;
    $params{parentname} = 'grp_pub';
    unless ($params{type_id}) {
	$params{cvname} = 'pubprop type' unless $params{cvname};
    }
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE grp_relationship element
# NOTE: this will create either a subject_id or object_id grp_relationship element
# but you have to attach this to the related grp elsewhere
# params
# doc - XML::DOM::Document required
# object _id - macro id for object grp or XML::DOM grp element
# subject_id - macro id for subject grp or XML::DOM grp element
# NOTE you can pass one or both of the above parameters with the following rules:
# if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other grp
# if both are passed then is_object, is_subject and any parameters to create a grp are ignored
# is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
# is_subject - boolean 't'          this flag indicates if the grp info provided should be 
#                                   added in as subject or object grp
# rtype_id - macro id for relationship type or XML::DOM cvterm element (Note: currently all is_relationship = '0'
# rtype - string for relationship type note: with this param  type will be assigned to relationship_type cv
# rank - integer optional with default = 0
# grp_id - macro_id for a grp or XML::DOM grp element required unless minimal grp bits provided
# uniquename - string required unless grp provided
# ftype_id - macro id for grp type or XML::DOM cvterm element required unless grp provided
# ftype - string for grp type (will be assigned to SO cv and can't be new as will be lookup)
sub create_ch_grp_relationship {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document
  
  ## grp_relationship element (will be returned)
  $params{relation} = $ldoc->createElement('grp_relationship');
  
  # deal with creating the grp bits if present and needed
  unless ($params{object_id} and $params{subject_id}) {
    # deal with creating the grp bits if present
    unless ($params{grp_id}) {
      # before creating grp with parameters
      print "WARNING - no required uniquename for grp\n" and return unless $params{uniquename};
 
      my %fparams = (doc => $ldoc,
		     uniquename => $params{uniquename},
		    );
      delete $params{uniquename};
      if ($params{ftype_id}) { 
	$fparams{type_id} = $params{ftype_id};
	delete $params{ftype_id}; 
      } elsif ($params{ftype}) {
	$fparams{type} = $params{ftype};
	delete $params{ftype};
      } else {
	print "WARNING -- you need to provide a grp type to make a grp!\n" and return;
	}
      $params{grp_id} = create_ch_grp(%fparams);
    } # now we have a grp element to associate as subject or object

    $params{thingy} = $params{grp_id};
    delete $params{grp_id};
  }
  # NOTE currently the cv for grp_relationship types is 'relationship type' 
  # uncomment and modify line below if this changes
  #$params{rtypecv} = 'relationship type' if ($params{rtype} and ! $params{rtypecv});
    return create_ch_relationship(%params);
}

# CREATE grp_relationshipprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from grp_relationshipprop type cv or XML::DOM cvterm element required 
#        Note: will default to making a grpprop from 'fr property type' cv unless cvname is provided
# cvname - string optional 
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_grp_relationshipprop {
    my %params = @_;
    $params{parentname} = 'grp_relationship';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE grp_relationship_pub element
# # Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of grp_relationship_pub
#      or just appending the pub element to grp_relationship_pub if that is passed
# params
# doc - XML::DOM::Document required
# grp_relationship_id -optional grp_relationship XML::DOM element or macro id
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_grp_relationship_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('grp_relationship_pub');

    if ($params{grp_relationship_id}) {
	$fp_el->appendChild(_build_element($ldoc,'grp_relationship_id',$params{grp_relationship_id}));
	delete $params{grp_relationship_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE grp_synonym element
# params
# doc - XML::DOM::Document required
# synonym_id - macro id for synonym? or XML::DOM synonym element required unless name and type provided
# name - string required unless synonym_id element provided
# type_id - macro id for synonym type or XML::DOM cvterm element
# type - string = name from the 'synonym type' cv
# pub_id - macro id for pub or a XML::DOM pub element required
# pub - a pub uniquename (i.e. FBrf)
# synonym_sgml - string optional but if not provided then synonym_sgml = name
#               - do not provide if a synonym element is provided
# is_current - optional string = 'f' or 't' default is 't' so don't provide this param
#               unless you know you want to change the value
# is_internal - optional string = 't' default is 'f' so don't provide this param
#               unless you know you want to change the value
sub create_ch_grp_synonym {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    ## grp_synonym element (will be returned)
    my $fs_el = $ldoc->createElement('grp_synonym');

    #create a synonym element and undefine name, type and synonym_sgml bits
    unless($params{synonym_id}) {
      # gather params for synonym
      print "WARNING - you haven't provided info to create a synonym element\n" and return unless $params{name};
      my %syn_params = (doc => $ldoc,
			name => $params{name},
		       );
      delete $params{name};
      if ($params{synonym_sgml}) {
	$syn_params{synonym_sgml} = $params{synonym_sgml};
	delete $params{synonym_sgml};
      }
       
      # check for type or type_id
      if ($params{type}) {
	$syn_params{type} = $params{type};
	delete $params{type};
      } elsif ($params{type_id}) {
	$syn_params{type_id} = $params{type_id};
	delete $params{type_id};
      } else {
	print "WARNING - you haven't provided a synonym type\n" and return;
      }
      my $syn_el = create_ch_synonym(%syn_params);
      $params{synonym_id} = $syn_el;	
    }

    # check for pub 
    if ($params{pub}) {
      $params{pub_id} = create_ch_pub(doc => $ldoc,
				      uniquename => $params{pub},
				     );
      delete $params{pub};
    }

     print "WARNING - you need to provide pub info for the synonym link\n" and return unless $params{pub_id};

    foreach my $e (keys %params) {
      next if ($e eq 'doc');
      $fs_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    
    return $fs_el;
}

# CREATE grpmember element
# params
# doc - XML::DOM::Document required
# grp_id - optional macro id for a grp or XML::DOM grp element for standalone grpmember
# type_id - macro id for grpmember type or XML::DOM cvterm element required 
# type - string from grpmember type cv
#        Note: will default to making a grpmember type from above cv unless cvname is provided
# cvname - string optional 
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value
# with_id - boolean optional if 1 then grpmember_id element is returned
# no_lookup - boolean option if 1 then default op="lookup" attribute will not be added to element
sub create_ch_grpmember {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $pdid_el = $ldoc->createElement('grpmember_id') if $params{with_id};

    ## grpmember element (will be returned)
    my $pd_el = $ldoc->createElement('grpmember');
    $pd_el->setAttribute('id',$params{macro_id}) if $params{macro_id};    

    # add an op="lookup" attribute unless no_lookup is specified
    unless ($params{no_lookup}) {
	$pd_el->setAttribute('op','lookup');
    }
    my $cv = 'grpmember type';
    if ($params{cvname}) {
	$cv = $params{cvname};
	delete $params{cvname};
    }

    # now deal make a cvterm element for type_id if string is provided
    if ($params{type}) {
      $params{type_id} = create_ch_cvterm(doc => $ldoc,
				   name => $params{type},
				   cv => $cv,
				  );
      delete $params{type};
    }

    $params{rank} = '0' unless $params{rank};

    foreach my $e (keys %params) {
	next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id' || $e eq 'no_lookup');
	$pd_el->appendChild(_build_element($ldoc,$e,$params{$e}));
  }

    return $pd_el;

}

# CREATE grpmember_cvterm element
# params
# doc - XML::DOM::Document required
# cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
# name - string required unless cvterm_id provided Note: a cvterm has a lookup by default cannot make a new cvterm 
#                                                        with this method
# cv_id - macro id for cv or XML::DOM cv element required if name and not cv
# cv - string for name of cv required if name and not cv_id
# is_not - optional boolean 't' or 'f' with default = 'f' so don't pass unless you know you want to change
sub create_ch_grpmember_cvterm {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fct_el = $ldoc->createElement('grpmember_cvterm');

    #create a cvterm element if necessary
    unless ($params{cvterm_id}) {
      print "WARNING -- you haven't provided params to make a cvterm, Sorry\n" and return unless ($params{name});
      my %cvtparams = (doc => $ldoc,
		       name => $params{name},
		      );
      delete $params{name};

      if ($params{cv_id}) {
	$cvtparams{cv_id} = $params{cv_id};
	delete $params{cv_id}; 
      } elsif ($params{cv}) {
	$cvtparams{cv} = $params{cv};
	delete $params{cv}; 
      } else {
	print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
      }
      $params{cvterm_id} = create_ch_cvterm(%cvtparams);      
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fct_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    return $fct_el;
}

# CREATE grpmemberprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id or XML::DOM cvterm element required 
#        Note: will default to making a grpmemberprop from 'property type' cv unless cvname is provided
# type - string from property type cv
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_grpmemberprop {
    my %params = @_;
    $params{parentname} = 'grpmember';
    unless ($params{type_id}) {
	$params{cvname} = 'property type' unless $params{cvname};
    }
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE grpmemberprop_pub element
# Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of grpmemberprop_pub
#      or just appending the pub element to grpmemberprop_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_grpmemberprop_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('grpmemberprop_pub');

    if ($params{grpmemberprop_id}) {
	$fp_el->appendChild(_build_element($ldoc,'grpmemberprop_id',$params{grpmemberprop_id}));
	delete $params{grpmemberprop_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE grpmember_pub element
# Note that this is just calling create_ch_pub setting with_id = 1 
#      and adding returned pub_id element as a child of grpmember_pub
#      or just appending the pub element to grpmember_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - macro_id for pub or XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
#        creating a new pub (i.e. not null value but not part of unique key
# type - string for type from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_grpmember_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('grpmember_pub');

    if ($params{grpmember_id}) {
	$fp_el->appendChild(_build_element($ldoc,'grpmember_id',$params{grpmember_id}));
	delete $params{grpmember_id};
    }
    
    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE humanhealth element
# params
# doc - XML::DOM::Document required
# uniquename - string required
# name - string optional
# organism_id - organism macro id or XML::DOM organism element required if no genus and species
# genus - string required if no organism
# species - string required if no organism
# dbxref_id - dbxref macro id or XML::DOM dbxref element required unless accession and db
# accession - string optional
# version - int optional
# db - string dbname optional
# is_obsolete - boolean 't' or 'f' default = 'f' optional
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value
# with_id - boolean optional if 1 then humanhealth_id element is returned
# no_lookup - boolean option if 1 then default op="lookup" attribute will not be added to element
sub create_ch_humanhealth {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fid_el = $ldoc->createElement('humanhealth_id') if $params{with_id};

    ## humanhealth element (will be returned)
    my $f_el = $ldoc->createElement('humanhealth');
    $f_el->setAttribute('id',$params{macro_id}) if $params{macro_id};    

    ## add an op="lookup" attribute unless no_lookup is specified
    #unless ($params{no_lookup}) {
#	$f_el->setAttribute('op','lookup');
#    }
	
    #create organism_id element if genus and species are provided
    unless ($params{organism_id}) {
	_add_orgid_param(\%params);
    }

    if ($params{accession}) {
      _add_dbxrefid_param(\%params);
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id' || $e eq 'no_lookup');
	$f_el->appendChild(_build_element($ldoc, $e,$params{$e}));
    }

    if ($fid_el) {
	$fid_el->appendChild($f_el);
	return $fid_el;
    }
    return $f_el;
}

# CREATE humanhealth_cvterm element
# params
# doc - XML::DOM::Document required
# cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
# name - string required unless cvterm_id provided Note: a cvterm has a lookup by default cannot make a new cvterm 
#                                                        with this method
# cv_id - macro id for cv or XML::DOM cv element required if name and not cv
# cv - string for name of cv required if name and not cv_id
# pub_id - macro id for pub or XML::DOM pub element required unless pub
# pub - string = pub uniquename Note: as pub has lookup option by default can't make a new pub using this param
# is_not - optional boolean 't' or 'f' with default = 'f' so don't pass unless you know you want to change
sub create_ch_humanhealth_cvterm {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fct_el = $ldoc->createElement('humanhealth_cvterm');

    #create a cvterm element if necessary
    unless ($params{cvterm_id}) {
      print "WARNING -- you haven't provided params to make a cvterm, Sorry\n" and return unless ($params{name});
      my %cvtparams = (doc => $ldoc,
		       name => $params{name},
		      );
      delete $params{name};

      if ($params{cv_id}) {
	$cvtparams{cv_id} = $params{cv_id};
	delete $params{cv_id}; 
      } elsif ($params{cv}) {
	$cvtparams{cv} = $params{cv};
	delete $params{cv}; 
      } else {
	print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
      }
      $params{cvterm_id} = create_ch_cvterm(%cvtparams);      
    }

    #create a pub element if necessary
    if ($params{pub}) {
      $params{pub_id} = create_ch_pub(doc => $ldoc,
				      uniquename => $params{pub},
				     );
      delete $params{pub};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fct_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    return $fct_el;
}

# CREATE humanhealth_cvtermprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from  humanhealth_cvtermprop type cv or XML::DOM cvterm element required 
#        Note: will default to making a humanhealthprop from 'humanhealth_cvtermprop type' cv unless 
#              cvname is provided
# cvname - string optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_humanhealth_cvtermprop {
    my %params = @_;
    $params{parentname} = 'humanhealth_cvterm';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}


# CREATE humanhealth_dbxref element
# params
# doc - XML::DOM::Document required
# humanhealth_id - macro humanhealth id or XML::DOM humanhealth element optionaal to create freestanding humanhealth_dbxref
# dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
# accession - string required unless dbxref_id provided
# db_id - macro db id or XML::DOM db element required unless dbxref_id provided
# db - string name of db
# version - string optional
# description - string optional
# macro_id - optional string to add as ID attribute value to humanhealth_dbxref
# is_current - string 't' or 'f' boolean default = 't' so don't pass unless
#              this shoud be changed
sub create_ch_humanhealth_dbxref {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fd_el = $ldoc->createElement('humanhealth_dbxref');
    $fd_el->setAttribute('id',$params{macro_id}) if $params{macro_id};

    if ($params{humanhealth_id}) {
	$fd_el->appendChild(_build_element($ldoc,'humanhealth_id',$params{humanhealth_id}));
	delete $params{humanhealth_id};
    }

    my $ic;
    if (exists($params{is_current})) { #assign value to a var and then remove from params
      if ($params{is_current}) {
	$ic = $params{is_current};
      } else {
	$ic = 'false';
      }
      delete $params{is_current};
    }
    
    # create a dbxref element if necessary
    unless ($params{dbxref_id}) {
      print "WARNING - missing required parameters, NO GO.\n" and return unless 
	($params{accession} and ($params{db_id} or $params{db}));
      if ($params{db}) {
	$params{db_id} = create_ch_db(doc => $ldoc,
				      name => $params{db},
				     );
	delete $params{db};
      }
      
	
      $params{dbxref_id} = create_ch_dbxref(%params);
    }

    $fd_el->appendChild(_build_element($ldoc,'dbxref_id',$params{dbxref_id})); #add dbxref element
    $fd_el->appendChild(_build_element($ldoc,'is_current',$ic)) if $ic;

    return $fd_el;
}

# CREATE humanhealth_dbxrefprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for property type or XML::DOM cvterm element required
# type -  string from  property type cv 
#        Note: will default to making a humanhealth_dbxrefprop from 'humanhealth_dbxrefprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'property type' if other than 'humanhealth_dbxrefprop type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_humanhealth_dbxrefprop {
    my %params = @_;
    $params{parentname} = 'humanhealth_dbxref';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}
# CREATE humanhealth_dbxrefprop_pub element
# Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of humanhealth_dbxrefprop
#      or just appending the pub element to humanhealth_dbxrefprop if that is passed
# params
# doc - XML::DOM::Document required
# humanhealth_dbxrefprop_id - XML::DOM humanhealth_dbxrefprop element
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_humanhealth_dbxrefprop_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('humanhealth_dbxrefprop_pub');

    if ($params{humanhealth_dbxrefprop_id}) {
        $fp_el->appendChild(_build_element($ldoc,'humanhealth_dbxrefprop_id',$params{humanhealth_dbxrefprop_id}));
        delete $params{humanhealth_dbxrefprop_id};
    }

    unless ($params{pub_id}) {
        $params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE humanhealth_feature element
# params
# doc - XML::DOM::Document required
# organism_id - macro id for organism or XML::DOM organism element
# genus - string
# species - string
# NOTE: you can use the generic paramaters in the following cases:
#       1.  you are only building either a humanhealth or feature element and not both
#       2.  or both humanhealth and feature have the same organism
#       otherwise use the prefixed parameters
# WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
# humanhealth_id - macro id for humanhealth or XML::DOM humanhealth element
# hh_uniquename - string humanhealth uniquename
# hh_organism_id - macro id for organism or XML::DOM organism element to link to humanhealth
# hh_genus
# hh_species
# feature_id - macro id for feature or XML::DOM feature element
# feat_uniquename
# feat_organism_id - macro id for organism or XML::DOM organism element to link to feature
# feat_genus
# feat_species
# feat_type_id
# feat_type
# pub_id - macroid or pub element
# pub - str pub uniquename
sub create_ch_humanhealth_feature {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  #have to think about the checks a bit

  my $lf_el = $ldoc->createElement('humanhealth_feature');

  # deal with feature bits
  if ($params{feat_uniquename}) {
    print "ERROR -- you don't have required parameters to make a feature, NO GO!\n" and return
      unless (($params{organism_id} or $params{feat_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{feat_genus} and $params{feat_species}))
	      and ($params{feat_type_id} or $params{feat_type}));

    unless ($params{feat_organism_id} or ($params{feat_genus} and $params{feat_species})) {
      $params{feat_organism_id} = $params{organism_id} if $params{organism_id};
      $params{feat_genus} = $params{genus} if $params{genus};
      $params{feat_species} = $params{species} if $params{species};
    }

    # gather all the feature parameters
    my %fparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /feat_(.+)/) {
	$fparams{$1} = $params{$p};
      }
    }

    $params{feature_id} = create_ch_feature(%fparams);
  }

  unless ($params{pub_id}) {
    print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
    $params{pub_id} = create_ch_pub(doc => $ldoc,
				    uniquename => $params{pub},
				   );
    delete $params{pub};
  }
  

  # likewise deal with the humanhealth bits
  if ($params{hh_uniquename}) {
    print "ERROR -- you don't have required parameters to make a humanhealth, NO GO!\n" and return
      unless (($params{organism_id} or $params{hh_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{hh_genus} and $params{hh_species})));

    unless ($params{hh_organism_id} or ($params{hh_genus} and $params{hh_species})) {
      $params{hh_organism_id} = $params{organism_id} if $params{organism_id};
      $params{hh_genus} = $params{genus} if $params{genus};
      $params{hh_species} = $params{species} if $params{species};
    }

    # gather all the humanhealth parameters
    my %lparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /hh_(.+)/) {
	$lparams{$1} = $params{$p};
      }
    }

    $params{humanhealth_id} = create_ch_humanhealth(%lparams);
  }

  # and then add the feature, humanhealth or both to the humanhealth_feature element
  $lf_el->appendChild(_build_element($ldoc,'humanhealth_id',$params{humanhealth_id})) if $params{humanhealth_id};
  $lf_el->appendChild(_build_element($ldoc,'feature_id',$params{feature_id})) if $params{feature_id};
  $lf_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id})) if $params{pub_id};

  return $lf_el;
}


# CREATE humanhealth_featureprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for humanhealth_featureprop type or XML::DOM cvterm element required
# type -  string from  humanhealth_featureprop type cv 
#        Note: will default to making a featureprop from 'humanhealth_featureprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'humanhealth_featureprop type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_humanhealth_featureprop {
    my %params = @_;
    $params{parentname} = 'humanhealth_feature';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE humanhealth_phenotype element
# params
# doc - XML::DOM::Document required
# humanhealth_id - optional macro id for humanhealth or XML::DOM humanhealth element
# phenotype_id - optional macro id for phenotype or XML::DOM phenotype element 
# uniquename - string required if no phenotype_id
# observable_id - macro id for observable or XML::DOM cvterm element optional
# attr_id - macro id for attr or XML::DOM cvterm element optional
# cvalue_id - macro id for cvalue or XML::DOM cvterm element optional
# assay_id - macro id for assay or XML::DOM cvterm element optional
# pub_id - macro id for pub or XML::DOM pub element
# pub - string uniquename for pub
# value - string optional
# macro_id - optional string to specify as ID attribute for phenotype

sub create_ch_humanhealth_phenotype {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  #have to think about the checks a bit

  my $lf_el = $ldoc->createElement('humanhealth_phenotype');

  # deal with phenotype bits
  unless ($params{phenotype_id}) {
    print "ERROR -- you don't have required parameters to make a phenotype, NO GO!\n" and return
      unless $params{uniquename};

    # gather all the feature parameters
    my %pparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      next if ($p eq 'humanhealth_id');
      $pparams{$p} = $params{$p};
    }

    $params{phenotype_id} = create_ch_phenotype(%pparams);
  }

  unless ($params{pub_id}) {
    print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
    $params{pub_id} = create_ch_pub(doc => $ldoc,
				    uniquename => $params{pub},
				   );
    delete $params{pub};
  }
  

  # and then add the feature, humanhealth or both to the humanhealth_feature element
  $lf_el->appendChild(_build_element($ldoc,'humanhealth_id',$params{humanhealth_id})) if $params{humanhealth_id};
  $lf_el->appendChild(_build_element($ldoc,'phenotype_id',$params{phenotype_id})) if $params{phenotype_id};
  $lf_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id})) if $params{pub_id};

  return $lf_el;
}

# CREATE humanhealth_phenotypeprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for humanhealth_phenotypeprop type or XML::DOM cvterm element required
# type -  string from  humanhealth_phenotypeprop type cv 
#        Note: will default to making a featureprop from 'humanhealth_phenotypeprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'humanhealth_phenotypeprop type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_humanhealth_phenotypeprop {
    my %params = @_;
    $params{parentname} = 'humanhealth_phenotype';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE humanhealth_pub element
# Note that this is just calling create_ch_pub setting with_id = 1 
#      and adding returned pub_id element as a child of humanhealth_pub
#      or just appending the pub element to humanhealth_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - macro_id for pub or XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
#        creating a new pub (i.e. not null value but not part of unique key
# type - string for type from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_humanhealth_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('humanhealth_pub');

    if ($params{humanhealth_id}) {
	$fp_el->appendChild(_build_element($ldoc,'humanhealth_id',$params{humanhealth_id}));
	delete $params{humanhealth_id};
    }
    
    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# just calling create_ch_prop with correct params to make desired prop
# note that there are cases here where the value is null
sub create_ch_humanhealth_pubprop {
    my %params = @_;
    $params{parentname} = 'humanhealth_pub';
    unless ($params{type_id}) {
	$params{cvname} = 'pubprop type' unless $params{cvname};
    }
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE humanhealth_relationship element
# NOTE: this will create either a subject_id or object_id humanhealth_relationship element
# but you have to attach this to the related humanhealth elsewhere
# params
# doc - XML::DOM::Document required
# object _id - macro id for object humanhealth or XML::DOM humanhealth element
# subject_id - macro id for subject humanhealth or XML::DOM humanhealth element
# NOTE you can pass one or both of the above parameters with the following rules:
# if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other humanhealth
# if both are passed then is_object, is_subject and any parameters to create a humanhealth are ignored
# is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
# is_subject - boolean 't'          this flag indicates if the humanhealth info provided should be 
#                                   added in as subject or object humanhealth
# rtype_id - macro id for relationship type or XML::DOM cvterm element (Note: currently all is_relationship = '0'
# rtype - string for relationship type note: with this param  type will be assigned to relationship_type cv unless rtypecv provided
# rtypecv - string for name of cv for type of relationship - defaults to 'relationship type'
# humanhealth_id - macro_id for a humanhealth or XML::DOM humanhealth element required unless minimal humanhealth bits provided
# uniquename - string required unless humanhealth provided
# organism_id - macro id for organism or XML::DOM organism element required unless humanhealth or (genus & species) provided
# genus - string required unless humanhealth or organism provided
# species - string required unless humanhealth or organism provided

sub create_ch_humanhealth_relationship {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document
  
  ## humanhealth_relationship element (will be returned)
  my $fr_el = $ldoc->createElement('humanhealth_relationship');
  $params{relation} = $fr_el;
  
  # deal with creating the humanhealth bits if present and needed
  unless ($params{object_id} and $params{subject_id}) {
    unless ($params{humanhealth_id}) {
      unless ($params{organism_id}) {
	_add_orgid_param(\%params);
      }
      
      # before creating humanhealth figure out which parameters
      print "WARNING - no required uniquename for humanhealth\n" and return unless $params{uniquename};
      my %fparams = (doc => $ldoc,
		     uniquename => $params{uniquename},
		     organism_id => $params{organism_id},
		    );
      delete $params{uniquename};
      delete $params{organism_id};
      $params{humanhealth_id} = create_ch_humanhealth(%fparams);
    } # now we have a humanhealth element to associate as subject or object
    
    $params{thingy} = $params{humanhealth_id};
    delete $params{humanhealth_id};
  } 
  #$params{rtypecv} = 'relationship type' if ($params{rtype} and ! $params{rtypecv});   
  return create_ch_relationship(%params);
}
  

# CREATE humanhealth_relationship_pub element
# # Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of humanhealth_relationship_pub
#      or just appending the pub element to humanhealth_relationship_pub if that is passed
# params
# doc - XML::DOM::Document required
# humanhealth_relationship_id -optional humanhealth_relationship XML::DOM element or macro id
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_humanhealth_relationship_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('humanhealth_relationship_pub');

    if ($params{humanhealth_relationship_id}) {
	$fp_el->appendChild(_build_element($ldoc,'humanhealth_relationship_id',$params{humanhealth_relationship_id}));
	delete $params{humanhealth_relationship_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}



# CREATE humanhealth_synonym element
# params
# doc - XML::DOM::Document required
# synonym_id - macro id for synonym? or XML::DOM synonym element required unless name and type provided
# name - string required unless synonym_id element provided
# type_id - macro id for synonym type or XML::DOM cvterm element
# type - string = name from the 'synonym type' cv
# pub_id - macro id for pub or a XML::DOM pub element required
# pub - a pub uniquename (i.e. FBrf)
# synonym_sgml - string optional but if not provided then synonym_sgml = name
#               - do not provide if a synonym element is provided
# is_current - optional string = 'f' or 't' default is 't' so don't provide this param
#               unless you know you want to change the value
# is_internal - optional string = 't' default is 'f' so don't provide this param
#               unless you know you want to change the value

sub create_ch_humanhealth_synonym {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $ls_el = $ldoc->createElement('humanhealth_synonym');

    # this is exactly the same thing as create_ch_feature_synonym with different name
    # so call that method then change parentage?
    my $el = create_ch_feature_synonym(%params);
    my @children = $el->getChildNodes;
    $ls_el->appendChild($_) for @children;

    return $ls_el;
}


# CREATE humanhealthprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id or XML::DOM cvterm element required 
#        Note: will default to making a humanhealthprop from 'property type' cv unless cvname is provided
# type - string from property type cv
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_humanhealthprop {
    my %params = @_;
    $params{parentname} = 'humanhealth';
    unless ($params{type_id}) {
	$params{cvname} = 'property type' unless $params{cvname};
    }
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}


# CREATE humanhealthprop_pub element
# Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of humanhealthprop_pub
#      or just appending the pub element to humanhealthprop_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_humanhealthprop_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('humanhealthprop_pub');

    if ($params{humanhealthprop_id}) {
	$fp_el->appendChild(_build_element($ldoc,'humanhealthprop_id',$params{humanhealthprop_id}));
	delete $params{humanhealthprop_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}
# CREATE interaction element
# params
# doc - XML::DOM::Document required
# uniquename - string required
# type_id - cvterm macro id or XML::DOM cvterm element required unless type
# type - string for a cvterm Note: will default to using 'PSI-MI' cv 
# unless a cvname is provided
# cvname - string optional to specify a cv other than 'PSI-MI' for the 
#         type_id do not use if providing a cvterm element or macro id to type_id
# description - string optional
# is_obsolete - boolean optional default = false
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value
# with_id - boolean optional if 1 then interaction_id element is returned
sub create_ch_interaction {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  print "WARNING -- missing required parameters\n" and return unless $params{uniquename}
    and ($params{type} or $params{type_id});

  my $iid_el = $ldoc->createElement('interaction_id') if $params{with_id};

  ## interaction element (will be returned)
  my $int_el = $ldoc->createElement('interaction');
  $int_el->setAttribute('id',$params{macro_id}) if $params{macro_id};

  # do we want to add a lookup?

  my $cv = 'PSI-MI';
  if ($params{cvname}) {
    $cv = $params{cvname};
    delete $params{cvname};
  }

  # now deal make a cvterm element for type_id if string is provided
  if ($params{type}) {
    $params{type_id} = create_ch_cvterm(doc => $ldoc,
					name => $params{type},
					cv => $cv,
				       );
    delete $params{type};
  }

  foreach my $e (keys %params) {
    next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id' || $e eq 'no_lookup');
    $int_el->appendChild(_build_element($ldoc, $e,$params{$e}));
  }

  if ($iid_el) {
    $iid_el->appendChild($int_el);
    return $iid_el;
  }
  return $int_el;
}



# CREATE interactionprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from interaction property type cv or XML::DOM cvterm element required
#        Note: will default to making a from 'interaction property type' cv unless cvname is provided
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_interactionprop {
    my %params = @_;
    $params{parentname} = 'interaction';
    unless ($params{type_id}) {
	$params{cvname} = 'interaction property type' unless $params{cvname};
    }
    my $ep_el = create_ch_prop(%params);
    return $ep_el;
}

# CREATE interactionprop_pub element
# Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of interactionprop_pub
#      or just appending the pub element to interactionprop_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_interactionprop_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return
unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('interactionprop_pub');

    if ($params{interactionyprop_id}) {
        $fp_el->appendChild(_build_element($ldoc,'interactionprop_id',$params{interactionprop_id}));
        delete $params{interactionprop_id};
    }

    unless ($params{pub_id}) {
        $params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE interaction_cvterm element
# params
# doc - XML::DOM::Document required
# cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
# name - string required unless cvterm_id provided Note: a cvterm has a lookup by default cannot make a new cvterm 
#                                                        with this method
# cv_id - macro id for cv or XML::DOM cv element required if name and not cv
# cv - string for name of cv required if name and not cv_id
sub create_ch_interaction_cvterm {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fct_el = $ldoc->createElement('interaction_cvterm');

    #create a cvterm element if necessary
    unless ($params{cvterm_id}) {
      print "WARNING -- you haven't provided params to make a cvterm, Sorry\n" and return unless ($params{name});
      my %cvtparams = (doc => $ldoc,
		       name => $params{name},
		      );
      delete $params{name};

      if ($params{cv_id}) {
	$cvtparams{cv_id} = $params{cv_id};
	delete $params{cv_id}; 
      } elsif ($params{cv}) {
	$cvtparams{cv} = $params{cv};
	delete $params{cv}; 
      } else {
	print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
      }
      $params{cvterm_id} = create_ch_cvterm(%cvtparams);      
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fct_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    return $fct_el;
}

# CREATE interaction_cvtermprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from  interaction_cvtermprop type cv or XML::DOM cvterm element required 
#        Note: will default to making a interactionprop from 'interaction_cvtermprop type' cv unless 
#              cvname is provided
# cvname - string optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_interaction_cvtermprop {
    my %params = @_;
    $params{parentname} = 'interaction_cvterm';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE interaction_cell_line element
# params
# doc - XML::DOM::Document required
# interaction_id - macro id for interaction or XML::DOM interaction element
# int_uniquename - string interaction uniquename
# int_type_id - macro id for interaction type or XML::DOM cvterm element 
# int_type - string for interaction type from interaction type cv
# cell_line_id - macro id for cell_line or XML::DOM cell_line element
# cell_uniquename - string cell_line uniquename
# organism_id - macro id for organism or XML::DOM organism element
# genus - string
# species - string
# NOTE: organism info is only required if you are building a cell_line element
# pub_id -  macro id for pub or XML::DOM pub element
# pub - pub uniquename
sub create_ch_interaction_cell_line {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  #have to think about the checks a bit

  my $lf_el = $ldoc->createElement('interaction_cell_line');

  # deal with interaction bits
  if ($params{int_uniquename}) {
    print "ERROR -- you don't have required parameters to make an interaction, NO GO!\n" and return
	      unless ($params{int_type_id} or $params{int_type});

    # gather all the interaction parameters
    my %fparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /int_(.+)/) {
	$fparams{$1} = $params{$p};
      }
    }
    $params{interaction_id} = create_ch_interaction(%fparams);
  }

  # likewise deal with the cell_line bits
  if ($params{cell_uniquename}) {
    print "ERROR -- you don't have required parameters to make a cell_line, NO GO!\n" and return
      unless ($params{organism_id} or ($params{genus} and $params{species}));

    # gather all the cell_line parameters
    my %lparams = (doc => $ldoc,
		   uniquename => $params{cell_uniquename},);

    if ($params{organism_id}) {
      $lparams{organism_id} = $params{organism_id};
    } else {
      $lparams{genus} = $params{genus};
      $lparams{species} = $params{species};
    }

    $params{cell_line_id} = create_ch_cell_line(%lparams);
  }

  #create a pub element if necessary
  unless ($params{pub_id}) {
    print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
    $params{pub_id} = create_ch_pub(doc => $ldoc,
				    uniquename => $params{pub},
				   );
  }


  # and then add the cell_line, interaction or both to the interaction_cell_line element
  $lf_el->appendChild(_build_element($ldoc,'interaction_id',$params{interaction_id})) if $params{interaction_id};
  $lf_el->appendChild(_build_element($ldoc,'cell_line_id',$params{cell_line_id})) if $params{cell_line_id};
  $lf_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));

  return $lf_el;
}

# CREATE interaction_expression element
# params
# doc - XML::DOM::Document required
# interaction_id - macro id for interaction or XML::DOM interaction element
# int_uniquename - string interaction uniquename
# int_type_id - macro id for interaction type or XML::DOM cvterm element 
# int_type - string for interaction type from interaction type cv
# expression_id - macro id for expression or XML::DOM expression element
# exp_uniquename - string expression uniquename
# pub_id -  macro id for pub or XML::DOM pub element
# pub - pub uniquename
sub create_ch_interaction_expression {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  #have to think about the checks a bit

  my $lf_el = $ldoc->createElement('interaction_expression');

  # deal with interaction bits
  if ($params{int_uniquename}) {
    print "ERROR -- you don't have required parameters to make an interaction, NO GO!\n" and return
	      unless ($params{int_type_id} or $params{int_type});

    # gather all the interaction parameters
    my %fparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /int_(.+)/) {
	$fparams{$1} = $params{$p};
      }
    }
    $params{interaction_id} = create_ch_interaction(%fparams);
  }

  # likewise deal with the expression bits
  unless ($params{expression_id}) {
    unless ($params{exp_uniquename}) {
      print "ERROR -- you don't have required parameters to make a expression, NO GO!\n";
      return;
    }
    $params{expression_id} = create_ch_expression(doc => $ldoc,
						  uniquename => $params{exp_uniquename},);
  }

  #create a pub element if necessary
  unless ($params{pub_id}) {
    print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
    $params{pub_id} = create_ch_pub(doc => $ldoc,
				    uniquename => $params{pub},
				   );
  }

  # and then add the expression, interaction or both to the interaction_expression element
  $lf_el->appendChild(_build_element($ldoc,'interaction_id',$params{interaction_id})) if $params{interaction_id};
  $lf_el->appendChild(_build_element($ldoc,'expression_id',$params{expression_id})) if $params{expression_id};
  $lf_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));

  return $lf_el;
}

# CREATE interaction_expressionprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for interaction_expressionprop type or XML::DOM cvterm element required
# type -  string from  interaction_cvtermprop type cv 
#        Note: will default to making a featureprop from 'expression_cvtermprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'interaction_expressionprop type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_interaction_expressionprop {
    my %params = @_;
    $params{parentname} = 'interaction_expression';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE interaction_pub element
# Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of interaction_pub
#      or just appending the pub element to interaction_pub if that is passed
# params
# doc - XML::DOM::Document required
# interaction_id - optional interaction element or macro_id to create freestanding ele
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_interaction_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('interaction_pub');

    if ($params{interaction_id}) {
	$fp_el->appendChild(_build_element($ldoc,'interaction_id',$params{interaction_id}));
	delete $params{interaction_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE library or library_id element
# params
# doc - XML::DOM::Document
# uniquename - string required
# type_id - macro id for library type or XML::DOM cvterm element required
# type - string to specify library type from 'FlyBase miscellaneous CV' cv
# organism_id - macro id for organism or XML::DOM organism element required if no genus and species
# genus - string required if no organism
# species - string required if no organism
# name - string optional
# macro_id - string to specify as id attribute
# with_id - boolean optional if 1 then db_id element is returned
sub create_ch_library {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document
    print "WARNING -- required params are missing - info for uniquename, type and organism are needed\n" and return
      unless ($params{uniquename} and ($params{type_id} or $params{type}) and ($params{organism_id} or ($params{genus} and $params{species})));

    my $libid_el = $ldoc->createElement('library_id') if $params{with_id};
    my $lib_el = $ldoc->createElement('library');
    $lib_el->setAttribute('id',$params{macro_id}) if $params{macro_id};  

    #create organism_id element if genus and species are provided
    unless ($params{organism_id}) {
	_add_orgid_param(\%params);
    }

    # create type_id element if type is specified
    if ($params{type}) {
      $params{type_id} = create_ch_cvterm(doc => $ldoc,
					  name => $params{type},
					  cv => 'FlyBase miscellaneous CV',
				  );
      delete $params{type};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id');
	$lib_el->appendChild(_build_element($ldoc, $e,$params{$e}));
    }


    if ($libid_el) {
	$libid_el->appendChild($lib_el);
	return $libid_el;
    }
    return $lib_el;
}

# CREATE library_cvterm element
# params
# doc - XML::DOM::Document required
# cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
# name - string required unless cvterm_id provided Note: a cvterm has a lookup by default cannot make a new cvterm 
#                                                        with this method
# cv_id - macro id for cv or XML::DOM cv element required if name and not cv
# cv - string for name of cv required if name and not cv_id
# pub_id - macro id for pub or XML::DOM pub element required unless pub
# pub - string = pub uniquename Note: as pub has lookup option by default can't make a new pub using this param
sub create_ch_library_cvterm {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fct_el = $ldoc->createElement('library_cvterm');

    #create a cvterm element if necessary
    unless ($params{cvterm_id}) {
      print "WARNING -- you haven't provided params to make a cvterm, Sorry\n" and return unless ($params{name});
      my %cvtparams = (doc => $ldoc,
		       name => $params{name},
		      );
      delete $params{name};

      if ($params{cv_id}) {
	$cvtparams{cv_id} = $params{cv_id};
	delete $params{cv_id}; 
      } elsif ($params{cv}) {
	$cvtparams{cv} = $params{cv};
	delete $params{cv}; 
      } else {
	print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
      }
      $params{cvterm_id} = create_ch_cvterm(%cvtparams);      
    }

    #create a pub element if necessary
    if ($params{pub}) {
      $params{pub_id} = create_ch_pub(doc => $ldoc,
				      uniquename => $params{pub},
				     );
      delete $params{pub};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fct_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    return $fct_el;
}

# CREATE library_cvtermprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for cvterm from library_cvtermprop type cv or XML::DOM cvterm element required
# type -  string from  library_cvtermprop type cv
#         Note: will default to making a cvterm from 'library_cvtermprop type' cv unless
#              cvname is provided
# cvname - string optional but see above for type and do not provide if passing type_id
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_library_cvtermprop {
    my %params = @_;
    $params{parentname} = 'library_cvterm';
    my $lp_el = create_ch_prop(%params);
    return $lp_el;
}

sub create_ch_library_dbxref {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fd_el = $ldoc->createElement('library_dbxref');

    if ($params{library_id}) {
	$fd_el->appendChild(_build_element($ldoc,'library_id',$params{library_id}));
	delete $params{library_id};
    }

    my $ic;
    if (exists($params{is_current})) { #assign value to a var and then remove from params
      if ($params{is_current}) {
	$ic = $params{is_current};
      } else {
	$ic = 'false';
      }
      delete $params{is_current};
    }
    
    # create a dbxref element if necessary
    unless ($params{dbxref_id}) {
      print "WARNING - missing required parameters, NO GO.\n" and return unless 
	($params{accession} and ($params{db_id} or $params{db}));
      if ($params{db}) {
	$params{db_id} = create_ch_db(doc => $ldoc,
				      name => $params{db},
				     );
	delete $params{db};
      }
      
	
      $params{dbxref_id} = create_ch_dbxref(%params);
    }

    $fd_el->appendChild(_build_element($ldoc,'dbxref_id',$params{dbxref_id})); #add dbxref element
    $fd_el->appendChild(_build_element($ldoc,'is_current',$ic)) if defined($ic);

    return $fd_el;
}

# CREATE library_dbxrefprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for property type or XML::DOM cvterm element required
# type -  string from  property type cv 
#        Note: will default to making a library_dbxrefprop from 'library_dbxrefprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'property type' if other than 'library_dbxrefprop type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_library_dbxrefprop {
    my %params = @_;
    $params{parentname} = 'library_dbxref';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE libraryprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from realtionship property type cv or XML::DOM cvterm element required 
#        Note: will default to making a featureprop from 'property type' cv unless cvname is provided
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_libraryprop {
    my %params = @_;
    $params{parentname} = 'library';
    unless ($params{type_id}) {
	$params{cvname} = 'property type' unless $params{cvname};
    }
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

sub create_ch_libraryprop_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('libraryprop_pub');

    if ($params{libraryprop_id}) {
	$fp_el->appendChild(_build_element($ldoc,'libraryprop_id',$params{libraryprop_id}));
	delete $params{libraryprop_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE library_expression element
# params
# doc - XML::DOM::Document required
# library_id - OPTIONAL macro library id or XML::DOM library element to create freestanding library_expression 
# expression_id - macro expression id or XML::DOM expression element - required unless uniquename provided
# uniquename - string required unless expression_id
# pub_id - macro pub id or XML::DOM pub element - required unless puname provided
# pub - string uniquename for pub (note will have lookup so can't create new pub here)
sub create_ch_library_expression {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fd_el = $ldoc->createElement('library_expression');

    # create a expression element if necessary
    unless ($params{expression_id}) {
      print "WARNING - missing required expression info, NO GO.\n" and return unless $params{uniquename};
      $params{expression_id} = create_ch_expression(doc => $ldoc,
						   uniquename => $params{uniquename},);
      delete $params{uniquename};
    }

    #create a pub element if necessary
    unless ($params{pub_id}) {
	print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
	$params{pub_id} = create_ch_pub(doc => $ldoc,
					uniquename => $params{pub},
				       );
	delete $params{pub};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fd_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }

    return $fd_el;
}

# CREATE library_expressionprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for library_expressionprop type or XML::DOM cvterm element required
# type -  string from  property type cv 
#        Note: will default to making a library_expressionprop from 'library_expressionprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'library_expression property type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_library_expressionprop {
    my %params = @_;
    $params{parentname} = 'library_expression';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE library_feature element
# params
# doc - XML::DOM::Document required
# organism_id - macro id for organism or XML::DOM organism element
# genus - string
# species - string
# NOTE: you can use the generic paramaters in the following cases:
#       1.  you are only building either a library or feature element and not both
#       2.  or both library and feature have the same organism
#       otherwise use the prefixed parameters
# WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
# library_id - macro id for library or XML::DOM library element
# lib_uniquename - string library uniquename
# lib_organism_id - macro id for organism or XML::DOM organism element to link to library
# lib_genus
# lib_species
# lib_type_id
# lib_type
# feature_id - macro id for feature or XML::DOM feature element
# feat_uniquename
# feat_organism_id - macro id for organism or XML::DOM organism element to link to feature
# feat_genus
# feat_species
# feat_type_id
# feat_type
sub create_ch_library_feature {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  #have to think about the checks a bit

  my $lf_el = $ldoc->createElement('library_feature');

  # deal with feature bits
  if ($params{feat_uniquename}) {
    print "ERROR -- you don't have required parameters to make a feature, NO GO!\n" and return
      unless (($params{organism_id} or $params{feat_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{feat_genus} and $params{feat_species}))
	      and ($params{feat_type_id} or $params{feat_type}));

    unless ($params{feat_organism_id} or ($params{feat_genus} and $params{feat_species})) {
      $params{feat_organism_id} = $params{organism_id} if $params{organism_id};
      $params{feat_genus} = $params{genus} if $params{genus};
      $params{feat_species} = $params{species} if $params{species};
    }

    # gather all the feature parameters
    my %fparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /feat_(.+)/) {
	$fparams{$1} = $params{$p};
      }
    }

    $params{feature_id} = create_ch_feature(%fparams);
  }

  # likewise deal with the library bits
  if ($params{lib_uniquename}) {
    print "ERROR -- you don't have required parameters to make a library, NO GO!\n" and return
      unless (($params{organism_id} or $params{lib_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{lib_genus} and $params{lib_species}))
	      and ($params{lib_type_id} or $params{lib_type}));

    unless ($params{lib_organism_id} or ($params{lib_genus} and $params{lib_species})) {
      $params{lib_organism_id} = $params{organism_id} if $params{organism_id};
      $params{lib_genus} = $params{genus} if $params{genus};
      $params{lib_species} = $params{species} if $params{species};
    }

    # gather all the library parameters
    my %lparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /lib_(.+)/) {
	$lparams{$1} = $params{$p};
      }
    }

    $params{library_id} = create_ch_library(%lparams);
  }

  # and then add the feature, library or both to the library_feature element
  $lf_el->appendChild(_build_element($ldoc,'library_id',$params{library_id})) if $params{library_id};
  $lf_el->appendChild(_build_element($ldoc,'feature_id',$params{feature_id})) if $params{feature_id};

  return $lf_el;
}

# CREATE library_featureprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for library_featureprop type or XML::DOM cvterm element required
# type -  string from  library_featureprop type cv 
#        Note: will default to making a featureprop from 'library_featureprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'library_featureprop type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_library_featureprop {
    my %params = @_;
    $params{parentname} = 'library_feature';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE library_grpmember element
# params 
# doc - XML::DOM::Document optional - required
#
# Here are parameters to make a grpmember element must either pass a grpmember_id or the other necessary bits
# grpmember_id - macro grpmember id or XML::DOM grpmember element
#
# Here are parameters to make an library element must either pass library_id or required bits
# library_id - macro library id or XML::DOM library element

sub create_ch_library_grpmember {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $af_el = $ldoc->createElement('library_grpmember');
     
    # or if a macro id or existing element have been provided
    if ($params{grpmember_id}) {
      $af_el->appendChild(_build_element($ldoc,'grpmember_id',$params{grpmember_id}));
      delete $params{grpmember_id};
    }
    if ($params{library_id}) {
      $af_el->appendChild(_build_element($ldoc,'library_id',$params{library_id}));
      delete $params{library_id};
    }

    return $af_el;
}

# CREATE library_humanhealth element
# params
# doc - XML::DOM::Document required
# organism_id - macro id for organism or XML::DOM organism element
# genus - string
# species - string
# NOTE: you can use the generic paramaters in the following cases:
#       1.  you are only building either a library or humanhealth element and not both
#       2.  or both library and humanhealth have the same organism
#       otherwise use the prefixed parameters
# WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
# library_id - macro id for library or XML::DOM library element
# lib_uniquename - string library uniquename
# lib_organism_id - macro id for organism or XML::DOM organism element to link to library
# lib_genus
# lib_species
# lib_type_id
# lib_type
# humanhealth_id - macro id for humanhealth or XML::DOM humanhealth element
# hh_uniquename
# hh_organism_id - macro id for organism or XML::DOM organism element to link to humanhealth
# hh_genus
# hh_species
# pub_id -  macro id for pub or XML::DOM pub element
# pub - pub uniquename
sub create_ch_library_humanhealth {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  #have to think about the checks a bit

  my $lf_el = $ldoc->createElement('library_humanhealth');

  # deal with humanhealth bits
  if ($params{hh_uniquename}) {
    print "ERROR -- you don't have required parameters to make a humanhealth, NO GO!\n" and return
      unless (($params{organism_id} or $params{hh_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{hh_genus} and $params{hh_species})));

    unless ($params{hh_organism_id} or ($params{hh_genus} and $params{hh_species})) {
      $params{hh_organism_id} = $params{organism_id} if $params{organism_id};
      $params{hh_genus} = $params{genus} if $params{genus};
      $params{hh_species} = $params{species} if $params{species};
    }

    # gather all the humanhealth parameters
    my %fparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /hh_(.+)/) {
	$fparams{$1} = $params{$p};
      }
    }

    $params{humanhealth_id} = create_ch_humanhealth(%fparams);
  }

  # likewise deal with the library bits
  if ($params{lib_uniquename}) {
    print "ERROR -- you don't have required parameters to make a library, NO GO!\n" and return
      unless (($params{organism_id} or $params{lib_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{lib_genus} and $params{lib_species}))
	      and ($params{lib_type_id} or $params{lib_type}));

    unless ($params{lib_organism_id} or ($params{lib_genus} and $params{lib_species})) {
      $params{lib_organism_id} = $params{organism_id} if $params{organism_id};
      $params{lib_genus} = $params{genus} if $params{genus};
      $params{lib_species} = $params{species} if $params{species};
    }

    # gather all the library parameters
    my %lparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /lib_(.+)/) {
	$lparams{$1} = $params{$p};
      }
    }

    $params{library_id} = create_ch_library(%lparams);
  }
  #create a pub element if necessary
  unless ($params{pub_id}) {
    print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
    $params{pub_id} = create_ch_pub(doc => $ldoc,
				    uniquename => $params{pub},
				   );
  }
  # and then add the humanhealth, library or both to the library_humanhealth element
  $lf_el->appendChild(_build_element($ldoc,'library_id',$params{library_id})) if $params{library_id};
  $lf_el->appendChild(_build_element($ldoc,'humanhealth_id',$params{humanhealth_id})) if $params{humanhealth_id};
  $lf_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));

  return $lf_el;
}

# CREATE library_humanhealthprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for library_humanhealthprop type or XML::DOM cvterm element required
# type -  string from  library_humanhealthprop type cv 
#        Note: will default to making a humanhealthprop from 'library_humanhealthprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'library_humanhealthprop type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_library_humanhealthprop {
    my %params = @_;
    $params{parentname} = 'library_humanhealth';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE library_interaction element
# params
# doc - XML::DOM::Document required
# library_id - macro id for library or XML::DOM library element
# lib_uniquename - string library uniquename
# organism_id - macro id for organism or XML::DOM organism element
# genus - string
# species - string
# NOTE: organism info is only required if you are building a library element
# lib_type_id - macro id for library type or XML::DOM cvterm element 
# lib_type - string for library type from FlyBase miscellaneous CV cv
# interaction_id - macro id for interaction or XML::DOM interaction element
# int_uniquename - string interaction uniquename
# int_type_id - macro id for interaction type or XML::DOM cvterm element
# int_type - string for interaction type from psi-mi ontology
# pub_id -  macro id for pub or XML::DOM pub element
# pub - pub uniquename
sub create_ch_library_interaction {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  #have to think about the checks a bit

  my $lf_el = $ldoc->createElement('library_interaction');

  # deal with interaction bits
  if ($params{int_uniquename}) {
    print "ERROR -- you don't have required parameters to make a interaction, NO GO!\n" and return
	      unless ($params{int_type_id} or $params{int_type});

    # gather all the interaction parameters
    my %fparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /int_(.+)/) {
	$fparams{$1} = $params{$p};
      }
    }

    $params{interaction_id} = create_ch_interaction(%fparams);
  }

  # likewise deal with the library bits
  if ($params{lib_uniquename}) {
    print "ERROR -- you don't have required parameters to make a library, NO GO!\n" and return
      unless (($params{organism_id} or ($params{genus} and $params{species}))
	      and ($params{lib_type_id} or $params{lib_type}));

    unless ($params{lib_organism_id} or ($params{lib_genus} and $params{lib_species})) {
      $params{lib_organism_id} = $params{organism_id} if $params{organism_id};
      $params{lib_genus} = $params{genus} if $params{genus};
      $params{lib_species} = $params{species} if $params{species};
    }

    # gather all the library parameters
    my %lparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /lib_(.+)/) {
	$lparams{$1} = $params{$p};
      }
    }

    $params{library_id} = create_ch_library(%lparams);
  }

  #create a pub element if necessary
  unless ($params{pub_id}) {
    print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
    $params{pub_id} = create_ch_pub(doc => $ldoc,
				    uniquename => $params{pub},
				   );
  }


  # and then add the interaction, library or both to the library_interaction element
  $lf_el->appendChild(_build_element($ldoc,'library_id',$params{library_id})) if $params{library_id};
  $lf_el->appendChild(_build_element($ldoc,'interaction_id',$params{interaction_id})) if $params{interaction_id};
  $lf_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));

  return $lf_el;
}

# CREATE library_pub element
# Note that this is just calling create_ch_pub  
#      and adding returned pub_id element as a child of library_pub
#      or just appending the pub element to library_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
#        creating a new pub (i.e. not null value but not part of unique key
# type - string for type from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_library_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('library_pub');
    if ($params{library_id}) {
	$fp_el->appendChild(_build_element($ldoc,'library_id',$params{library_id}));
	delete $params{library_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE library_relationship element
# NOTE: this will create either a subject_id or object_id library_relationship element
# but you have to attach this to the related library elsewhere
# params
# doc - XML::DOM::Document required
# object _id - macro id for object library or XML::DOM library element
# subject_id - macro id for subject library or XML::DOM library element
# NOTE you can pass one or both of the above parameters with the following rules:
# if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other library
# if both are passed then is_object, is_subject and any parameters to create a library are ignored
# is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
# is_subject - boolean 't'          this flag indicates if the library info provided should be 
#                                   added in as subject or object library
# rtype_id - macro id for relationship type or XML::DOM cvterm element (Note: currently all is_relationship = '0'
# rtype - string for relationship type note: with this param  type will be assigned to relationship_type cv
# library_id - macro_id for a library or XML::DOM library element required unless minimal library bits provided
# uniquename - string required unless library provided
# organism_id - macro id for organism or XML::DOM organism element required unless library or (genus & species) provided
# genus - string required unless library or organism provided
# species - string required unless library or organism provided
# ftype_id - macro id for library type or XML::DOM cvterm element required unless library provided
# ftype - string for library type
sub create_ch_library_relationship {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document
  
  ## library_relationship element (will be returned)
  $params{relation} = $ldoc->createElement('library_relationship');
  
  # deal with creating the library bits if present
  unless ($params{object_id} and $params{subject_id}) {
    unless ($params{library_id}) {
      unless ($params{organism_id}) {
	_add_orgid_param(\%params);
      }
	
      # before creating strain library which parameters
      print "WARNING - no required uniquename for library\n" and return unless $params{uniquename};
      my %fparams = (doc => $ldoc,
		     uniquename => $params{uniquename},
		     organism_id => $params{organism_id},
		    );
      delete $params{uniquename};
      delete $params{organism_id};
      if ($params{stype_id}) { 
	$fparams{type_id} = $params{ftype_id};
	delete $params{ftype_id};
      } elsif ($params{ftype}) {
	$fparams{type} = $params{ftype};
	delete $params{ftype};
      } else {
	print "WARNING -- you need to provide a library type to make a library!\n" and return;
      }
      $params{library_id} = create_ch_library(%fparams);
    } # now we have a strain element to associate as subject or object
    
    $params{thingy} = $params{library_id};
    delete $params{library_id};
  } 
  # NOTE currently the cv for library_relationship types is 'relationship type' 
  # uncomment and modify line below if this changes
  #$params{rtypecv} = 'relationship type' if ($params{rtype} and ! $params{rtypecv});
  return create_ch_relationship(%params);
}
  

# CREATE library_relationship_pub element
# # Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of library_relationship_pub
#      or just appending the pub element to library_relationship_pub if that is passed
# params
# doc - XML::DOM::Document required
# library_relationship_id -optional library_relationship XML::DOM element or macro id
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_library_relationship_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('library_relationship_pub');

    if ($params{library_relationship_id}) {
	$fp_el->appendChild(_build_element($ldoc,'library_relationship_id',$params{library_relationship_id}));
	delete $params{library_relationship_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE library_strain element
# params
# doc - XML::DOM::Document required
# organism_id - macro id for organism or XML::DOM organism element
# genus - string
# species - string
# NOTE: you can use the generic paramaters in the following cases:
#       1.  you are only building either a library or strain element and not both
#       2.  or both library and strain have the same organism
#       otherwise use the prefixed parameters
# WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
# library_id - macro id for library or XML::DOM library element
# lib_uniquename - string library uniquename
# lib_organism_id - macro id for organism or XML::DOM organism element to link to library
# lib_genus
# lib_species
# lib_type_id
# lib_type
# strain_id - macro id for strain or XML::DOM strain element
# str_uniquename
# str_organism_id - macro id for organism or XML::DOM organism element to link to strain
# str_genus
# str_species
# pub_id -  macro id for pub or XML::DOM pub element
# pub - pub uniquename
sub create_ch_library_strain {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  #have to think about the checks a bit

  my $lf_el = $ldoc->createElement('library_strain');

  # deal with strain bits
  if ($params{str_uniquename}) {
    print "ERROR -- you don't have required parameters to make a strain, NO GO!\n" and return
      unless (($params{organism_id} or $params{str_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{str_genus} and $params{str_species})));

    unless ($params{str_organism_id} or ($params{str_genus} and $params{str_species})) {
      $params{str_organism_id} = $params{organism_id} if $params{organism_id};
      $params{str_genus} = $params{genus} if $params{genus};
      $params{str_species} = $params{species} if $params{species};
    }

    # gather all the strain parameters
    my %fparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /str_(.+)/) {
	$fparams{$1} = $params{$p};
      }
    }

    $params{strain_id} = create_ch_strain(%fparams);
  }

  # likewise deal with the library bits
  if ($params{lib_uniquename}) {
    print "ERROR -- you don't have required parameters to make a library, NO GO!\n" and return
      unless (($params{organism_id} or $params{lib_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{lib_genus} and $params{lib_species}))
	      and ($params{lib_type_id} or $params{lib_type}));

    unless ($params{lib_organism_id} or ($params{lib_genus} and $params{lib_species})) {
      $params{lib_organism_id} = $params{organism_id} if $params{organism_id};
      $params{lib_genus} = $params{genus} if $params{genus};
      $params{lib_species} = $params{species} if $params{species};
    }

    # gather all the library parameters
    my %lparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /lib_(.+)/) {
	$lparams{$1} = $params{$p};
      }
    }

    $params{library_id} = create_ch_library(%lparams);
  }
  #create a pub element if necessary
  unless ($params{pub_id}) {
    print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
    $params{pub_id} = create_ch_pub(doc => $ldoc,
				    uniquename => $params{pub},
				   );
  }
  # and then add the strain, library or both to the library_strain element
  $lf_el->appendChild(_build_element($ldoc,'library_id',$params{library_id})) if $params{library_id};
  $lf_el->appendChild(_build_element($ldoc,'strain_id',$params{strain_id})) if $params{strain_id};
  $lf_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));

  return $lf_el;
}

# CREATE library_strainprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for library_strainprop type or XML::DOM cvterm element required
# type -  string from  library_strainprop type cv 
#        Note: will default to making a strainprop from 'library_strainprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'library_strainprop type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_library_strainprop {
    my %params = @_;
    $params{parentname} = 'library_strain';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}


# CREATE library_synonym element
# params
# doc - XML::DOM::Document required
# synonym_id - XML::DOM synonym element required unless name and type provided
# name - string required unless synonym_id element provided
# type_id - macro id for synonym type or XML::DOM cvterm element
#                 required unless a synonym element provided
# type - string = name from the 'synonym type' cv
# pub_id macro id for a pub or a XML::DOM pub element required
# pub - a pub uniquename (i.e. FBrf)
# synonym_sgml - string optional but if not provided then synonym_sgml = name
#               - do not provide if a synonym element is provided
# is_current - optional string = 'f' or 't' default is 't' so don't provide this param
#               unless you know you want to change the value
# is_internal - optional string = 't' default is 'f' so don't provide this param
#               unless you know you want to change the value
sub create_ch_library_synonym {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $ls_el = $ldoc->createElement('library_synonym');

    # this is exactly the same thing as create_ch_feature_synonym with different name
    # so call that method then change parentage?
    my $el = create_ch_feature_synonym(%params);
    my @children = $el->getChildNodes;
    $ls_el->appendChild($_) for @children;

    return $ls_el;
}

# CREATE organism or organism_id element
# params
# doc - XML::DOM::Document
# genus - string required
# species - string required
# abbreviation - string optional
# common_name - string optional
# comment - string optional
# macro_id - string to specify with id attribute
# with_id - boolean optional if 1 then organism_id element is returned
sub create_ch_organism {
    my %params = @_;
    $params{elname} = 'organism';
    $params{required} = ['genus','species'];
    my $eel = _create_simple_element(%params);
    return $eel;
}

# CREATE organism_cvterm element
# params
# doc - XML::DOM::Document required
# cvterm_id - XML::DOM cvterm element unless other cvterm bits are provided
# name - string required unless cvterm_id provided
# cv_id -  macro id for a cv or XML::DOM cv element required unless cvterm_id or cv provided
# cv - string = cvname required unless cvterm_id or cv_id provided
# organism_id - macro id for organism or organism element if you want to create a freestanding organism_cvterm
sub create_ch_organism_cvterm {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  my $ect_el = $ldoc->createElement('organism_cvterm');

  # create a cvterm element if necessary
  unless ($params{cvterm_id}) {
    print "WARNING -- you haven't provided params to make a cvterm, Sorry\n" and return unless ($params{name});
    my %cvtparams = (doc => $ldoc,
		     name => $params{name},
		    );
    delete $params{name};
    
    if ($params{cv_id}) {
      $cvtparams{cv_id} = $params{cv_id};
      delete $params{cv_id}; 
    } elsif ($params{cv}) {
      $cvtparams{cv} = $params{cv};
      delete $params{cv}; 
    } else {
      print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
    }
    $params{cvterm_id} = create_ch_cvterm(%cvtparams);      
  }

  #create a pub element if necessary
    if ($params{pub}) {
      $params{pub_id} = create_ch_pub(doc => $ldoc,
				      uniquename => $params{pub},
				     );
      delete $params{pub};
    }
  
  #now set required rank to 0 if not provided
  $params{rank} = '0' unless $params{rank};

  foreach my $e (keys %params) {
    next if ($e eq 'doc');
    $ect_el->appendChild(_build_element($ldoc,$e,$params{$e}));
  }
  return $ect_el;
}

# CREATE organism_cvtermprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for cvterm from organism_cvtprop type cv or XML::DOM cvterm element required
# type -  string from  organism_cvtermprop type cv
#         Note: will default to making a cvterm from 'organism_cvtermprop type' cv unless
#              cvname is provided
# cvname - string optional but see above for type and do not provide if passing type_id
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_organism_cvtermprop {
    my %params = @_;
    $params{parentname} = 'organism_cvterm';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE organism_dbxref element
# params
# doc - XML::DOM::Document required
# organism_id - macro organism id or XML::DOM organism element optionaal to create freestanding organism_dbxref
# dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
# accession - string required unless dbxref_id provided
# db_id - macro db id or XML::DOM db element required unless dbxref_id provided
# db - string name of db
# version - string optional
# description - string optional
# is_current - string 't' or 'f' boolean default = 't' so don't pass unless
#              this shoud be changed
sub create_ch_organism_dbxref {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fd_el = $ldoc->createElement('organism_dbxref');

    if ($params{organism_id}) {
	$fd_el->appendChild(_build_element($ldoc,'organism_id',$params{organism_id}));
	delete $params{organism_id};
    }

    my $ic;
    if (exists($params{is_current})) { #assign value to a var and then remove from params
      if ($params{is_current}) {
	$ic = $params{is_current};
      } else {
	$ic = 'false';
      }
      delete $params{is_current};
    }
    
    # create a dbxref element if necessary
    unless ($params{dbxref_id}) {
      print "WARNING - missing required parameters, NO GO.\n" and return unless 
	($params{accession} and ($params{db_id} or $params{db}));
      if ($params{db}) {
	$params{db_id} = create_ch_db(doc => $ldoc,
				      name => $params{db},
				     );
	delete $params{db};
      }
      
	
      $params{dbxref_id} = create_ch_dbxref(%params);
    }

    $fd_el->appendChild(_build_element($ldoc,'dbxref_id',$params{dbxref_id})); #add dbxref element
    $fd_el->appendChild(_build_element($ldoc,'is_current',$ic)) if $ic;

    return $fd_el;
}


# CREATE organism_pub element
# Note that this is just calling create_ch_pub setting with_id = 1
#      and adding returned pub_id element as a child of organism_pub
#      or just appending the pub element to organism_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
#            creating a new pub (i.e. not null value but not part of unique key
# type - string from pub type cv same requirement as type_id
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_organism_pub {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document
  
  my $ep_el = $ldoc->createElement('organism_pub');
  
    if ($params{organism_id}) {
	$ep_el->appendChild(_build_element($ldoc,'organism_id',$params{organism_id}));
	delete $params{organism_id};
    }

  unless ($params{pub_id}) {
    $params{pub_id} = create_ch_pub(%params); #will return a pub element
  }
  $ep_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));

  return $ep_el;
}

# CREATE organismprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from realtionship property type cv or XML::DOM cvterm element required
#        Note: will default to making a organismprop from 'property type' cv unless cvname is provided
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_organismprop {
    my %params = @_;
    $params{parentname} = 'organism';
    unless ($params{type_id}) {
	$params{cvname} = 'organismprop type' unless $params{cvname};
    }
    my $ep_el = create_ch_prop(%params);
    return $ep_el;
}

# CREATE organismprop_pub element
# # Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of organismprop_pub
#      or just appending the pub element to organismprop_pub if that is passed
# params
# doc - XML::DOM::Document required
# organismprop_id -optional organismprop XML::DOM element or macro id
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_organismprop_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('organismprop_pub');

    if ($params{organismprop_id}) {
	$fp_el->appendChild(_build_element($ldoc,'organismprop_id',$params{organismprop_id}));
	delete $params{organismprop_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE organism_grpmember element
# params 
# doc - XML::DOM::Document optional - required
#
# Here are parameters to make a grpmember element must either pass a grpmember_id or the other necessary bits
# grpmember_id - macro grpmember id or XML::DOM grpmember element
#
# Here are parameters to make an organism element must either pass organism_id or required bits
# organism_id - macro organism id or XML::DOM organism element

sub create_ch_organism_grpmember {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $af_el = $ldoc->createElement('organism_grpmember');
     
    # or if a macro id or existing element have been provided
    if ($params{grpmember_id}) {
      $af_el->appendChild(_build_element($ldoc,'grpmember_id',$params{grpmember_id}));
      delete $params{grpmember_id};
    }

   # and here we are dealing with organism info if provided
    if ($params{organism_id}) {
      $af_el->appendChild(_build_element($ldoc,'organism_id',$params{organism_id}));
      delete $params{organism_id};
    }

    return $af_el;
}


# CREATE organism_library element
# params 
# doc - XML::DOM::Document optional - required
#
# Here are parameters to make a grpmember element must either pass a grpmember_id or the other necessary bits
# grpmember_id - macro grpmember id or XML::DOM grpmember element
#
# Here are parameters to make an organism element must either pass organism_id or required bits
# organism_id - macro organism id or XML::DOM organism element

sub create_ch_organism_library {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $ol_el = $ldoc->createElement('organism_library');
     
    # or if a macro id or existing element have been provided
    if ($params{library_id}) {
      $ol_el->appendChild(_build_element($ldoc,'library_id',$params{library_id}));
      delete $params{library_id};
    }

   # and here we are dealing with organism info if provided
    if ($params{organism_id}) {
      $ol_el->appendChild(_build_element($ldoc,'organism_id',$params{organism_id}));
      delete $params{organism_id};
    }

    return $ol_el;
}


# CREATE phendesc element
# params
# doc - XML::DOM::Document required
# genotype_id - macro id for genotype or XML::DOM genotype element
# genotype - string genotype uniquename 
# environment_id - macro id for environment or XML::DOM environment element
# environment - string environment uniquename
# description - string optional but can't be null if creating new
# type_id - macro id for phendesc type or XML::DOM cvterm element
# type - string for cvterm name from phendesc type CV
# pub_id - macro id for pub or XML::DOM pub element
# pub - string pub uniquename
sub create_ch_phendesc {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document
    print "ERROR -- required parameters are missing: description, type or pub info - NO GO!\n" and return
      unless (($params{type_id} or $params{type}) and ($params{pub_id} or $params{pub}));
    print "ERROR -- you need either genotype info or environment info or both - NO GO!\n" and return
      unless ($params{genotype_id} or $params{genotype} or $params{environment_id} or $params{environment});


    my $pd_el = $ldoc->createElement('phendesc');  

    if ($params{genotype}) {
      $params{genotype_id} = create_ch_genotype(doc => $ldoc, uniquename => $params{genotype});
      delete $params{genotype};
    }

    if ($params{environment}) {
      $params{environment_id} = create_ch_environment(doc => $ldoc, uniquename => $params{environment});
      delete $params{environment};
    }

    if ($params{type}) {
      $params{type_id} = create_ch_cvterm(doc => $ldoc, name => $params{type}, cv => 'phendesc type');
      delete $params{type};
    }

    if ($params{pub}) {
      $params{pub_id} = create_ch_pub(doc => $ldoc, uniquename => $params{pub},);
      delete $params{pub};
    }

    foreach my $e (keys %params) {
      next if ($e eq 'doc');
	$pd_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }

    return $pd_el;
}

# CREATE phenotype or phenotype_id element
# params
# doc - XML::DOM::Document required
# uniquename - string required
# observable_id - macro id for observable or XML::DOM cvterm element optional
# attr_id - macro id for attr or XML::DOM cvterm element optional
# cvalue_id - macro id for cvalue or XML::DOM cvterm element optional
# assay_id - macro id for assay or XML::DOM cvterm element optional
# value - string optional
# macro_id - optional string to specify as ID attribute for genotype
# with_id - boolean optional if 1 then genotype_id element is returned
sub create_ch_phenotype {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document
    print "WARNING -- no uniquename specified\n" and return unless $params{uniquename};

    my $phid_el = $ldoc->createElement('phenotype_id') if $params{with_id};
    my $ph_el = $ldoc->createElement('phenotype');
    $ph_el->setAttribute('id',$params{macro_id}) if $params{macro_id};

    foreach my $e (keys %params) {
      next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id');
	$ph_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }

    if ($phid_el) {
	$phid_el->appendChild($ph_el);
	return $phid_el;
    }
    return $ph_el;
}

# CREATE phenotype_comparison element
# params
# doc - XML::DOM::Document required
# organism_id - macro id for an organism or XML::DOM organism element required unless genus and species
# genus - string required if no organism_id
# species - string required if no organism_id
# genotype1_id - macro id for a genotype or XML::DOM genotype element required unless genotype1
# genotype1 - string genotype uniquename required unless genotype1_id
# environment1_id - macro id for a environment or XML::DOM environment element required unless environment1
# environment1 - string environment uniquename required unless environment1_id
# genotype2_id - macro id for a genotype or XML::DOM genotype element required unless genotype2
# genotype2 - string genotype uniquename required unless genotype2_id
# environment2_id - macro id for a environment or XML::DOM environment element required unless environment2
# environment2 - string environment uniquename required unless environment2_id
# phenotype1_id - macro id for phenotype or XML::DOM phenotype element required unless phenotype1
# phenotype1 - string phenotype uniquename required unless phenotype1_id
# phenotype2_id - macro id for phenotype or XML::DOM phenotype element optional
# phenotype2 - string phenotype uniquename optional
# pub_id macro id for a pub or a XML::DOM pub element required unless pub
# pub - a pub uniquename (i.e. FBrf) required unless pub_id
sub create_ch_phenotype_comparison {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document
  print "ERROR -- required parameters are missing: organism, genotype1, environment1, genotype2, environment2, phenotype1 or pub info - NO GO!\n" 
    and return
    unless (($params{genotype1_id} or $params{genotype1}) and ($params{genotype2_id} or $params{genotype2})
	    and ($params{environment1_id} or $params{environment1}) and ($params{environment2_id} or $params{environment2})
	    and ($params{phenotype1_id} or $params{phenotype1}) and ($params{pub_id} or $params{pub}) 
	    and ($params{organism_id} or ($params{genus} and $params{species})));

  my $pc_el = $ldoc->createElement('phenotype_comparison');

  #create organism_id element if genus and species are provided
  unless ($params{organism_id}) {
 	_add_orgid_param(\%params);
  }

  if ($params{genotype1}) {
    $params{genotype1_id} = create_ch_genotype(doc => $ldoc, uniquename => $params{genotype1});
    delete $params{genotype1};
  }

  if ($params{environment1}) {
    $params{environment1_id} = create_ch_environment(doc => $ldoc, uniquename => $params{environment1});
    delete $params{environment1};
  }

  if ($params{genotype2}) {
    $params{genotype2_id} = create_ch_genotype(doc => $ldoc, uniquename => $params{genotype2});
    delete $params{genotype2};
  }

  if ($params{environment2}) {
    $params{environment2_id} = create_ch_environment(doc => $ldoc, uniquename => $params{environment2});
    delete $params{environment2};
  }

  if ($params{phenotype1}) {
    $params{phenotype1_id} = create_ch_phenotype(doc => $ldoc, uniquename => $params{phenotype1});
    delete $params{phenotype1};
  }

  if ($params{phenotype2}) {
    $params{phenotype2_id} = create_ch_phenotype(doc => $ldoc, uniquename => $params{phenotype2});
    delete $params{phenotype2};
  }

  if ($params{pub}) {
    $params{pub_id} = create_ch_pub(doc => $ldoc, uniquename => $params{pub},);
    delete $params{pub};
  }

  foreach my $e (keys %params) {
    next if ($e eq 'doc');
    $pc_el->appendChild(_build_element($ldoc,$e,$params{$e}));
  }

  return $pc_el;
}
*create_ch_ph_comp = \&create_ch_phenotype_comparison;

# CREATE phenotype_comparison_cvterm element
# params
# doc - XML::DOM::Document optional - required
# phenotype_comparison_id - optional macro id for phenotype_comparison or phenotype_comparison XML::DOM element
# NOTE: to make a standalone element
# cvterm_id -  macro id for cvterm of XML::DOM cvterm element
# name - cvterm name
# cv_id - macro id for a CV or XML::DOM cv element
# cv - name of a cv
# is_obsolete - optional param for cvterm
# NOTE: you need to pass cvterm bits if attaching to existing phenotype element or 
#       creating a freestanding phenotype_cvterm
# rank - optional with default = 0 so only pass if you want a different rank
sub create_ch_phenotype_comparison_cvterm {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  my $ect_el = $ldoc->createElement('phenotype_comparison_cvterm');

  # create a cvterm element if necessary
  if ($params{name}) {
    print "ERROR: You don't have all the parameters required to make a cvterm, NO GO!\n" and return
      unless ($params{cv_id} or $params{cv});
    my %cvtparams = (doc => $ldoc,
		     name => $params{name},
		    );
    delete $params{name};

    if ($params{cv_id}) {
      $cvtparams{cv_id} = $params{cv_id};
      delete $params{cv_id}; 
    } elsif ($params{cv}) {
      $cvtparams{cv} = $params{cv};
      delete $params{cv}; 
    } else {
      print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
    }

    if ($params{is_obsolete}) {
      $cvtparams{is_obsolete} = $params{is_obsolete};
      delete $params{is_obsolete};
    }
    $params{cvterm_id} = create_ch_cvterm(%cvtparams);
  }

  # now see which elements to attach to phenotype_cvterm
  $ect_el->appendChild(_build_element($ldoc,'phenotype_comparison_id',$params{phenotype_comparison_id})) if $params{phenotype_comparison_id};
  $ect_el->appendChild(_build_element($ldoc,'cvterm_id',$params{cvterm_id})) if $params{cvterm_id};

  #add the required rank element 
  my $rank = '0';
  $rank = $params{rank} if $params{rank};
  $ect_el->appendChild(_build_element($ldoc,'rank',$rank));

  return $ect_el;
}
*create_ch_ph_comp_cvt = \&create_ch_phenotype_comparison_cvterm;

# CREATE phenotype_cvterm element
# params
# doc - XML::DOM::Document optional - required
# phenotype_id - macro id for phenotype or XML::DOM phenotype element
# uniquename - phenotype uniquename
# NOTE: you need to pass phenotype bits if attaching to existing cvterm element or 
#       creating a freestanding phenotype_cvterm
# cvterm_id -  macro id for cvterm of XML::DOM cvterm element
# name - cvterm name
# cv_id - macro id for a CV or XML::DOM cv element
# cv - name of a cv
# is_obsolete - optional param for cvterm
# rank - optional with default = 0 so only pass if you want a different rank
# NOTE: you need to pass cvterm bits if attaching to existing phenotype element or 
#       creating a freestanding phenotype_cvterm
sub create_ch_phenotype_cvterm {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  my $ect_el = $ldoc->createElement('phenotype_cvterm');

  # create an phenotype if necessary
  if ($params{uniquename}) {
    $params{phenotype_id} = create_ch_phenotype(doc => $ldoc,
						uniquename => $params{uniquename},
					       );
    delete $params{uniquename};
  }

  # create a cvterm element if necessary
  if ($params{name}) {
    print "ERROR: You don't have all the parameters required to make a cvterm, NO GO!\n" and return
      unless ($params{cv_id} or $params{cv});
    my %cvtparams = (doc => $ldoc,
		     name => $params{name},
		    );
    delete $params{name};

    if ($params{cv_id}) {
      $cvtparams{cv_id} = $params{cv_id};
      delete $params{cv_id}; 
    } elsif ($params{cv}) {
      $cvtparams{cv} = $params{cv};
      delete $params{cv}; 
    } else {
      print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
    }

    if ($params{is_obsolete}) {
      $cvtparams{is_obsolete} = $params{is_obsolete};
      delete $params{is_obsolete};
    }
    $params{cvterm_id} = create_ch_cvterm(%cvtparams);      
  }

  # now see which elements to attach to phenotype_cvterm
  $ect_el->appendChild(_build_element($ldoc,'phenotype_id',$params{phenotype_id})) if $params{phenotype_id};
  $ect_el->appendChild(_build_element($ldoc,'cvterm_id',$params{cvterm_id})) if $params{cvterm_id};

  #add the required rank element 
  my $rank = '0';
  $rank = $params{rank} if $params{rank};
  $ect_el->appendChild(_build_element($ldoc,'rank',$rank));
  
  return $ect_el;
}
*create_ch_ph_cvt = \&create_ch_phenotype_cvterm;

# CREATE phenstatement element
# params
# doc - XML::DOM::Document required
# genotype_id - macro id for a genotype or XML::DOM genotype element required unless genotype
# genotype - string genotype uniquename required unless genotype_id
# environment_id - macro id for a environment or XML::DOM environment element required unless environment
# environment - string environment uniquename required unless environment_id
# phenotype_id - macro id for phenotype or XML::DOM phenotype element required unless phenotype
# phenotype - string phenotype uniquename required unless phenotype_id
# type_id - macro id for a phenstatement type or XML::DOM cvterm element
# pub_id - macro id for a pub or a XML::DOM pub element required unless pub
# pub - a pub uniquename (i.e. FBrf) required unless pub_id
sub create_ch_phenstatement {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document
  print "ERROR -- required parameters are missing: genotype, environment, phenotype, type and pub all required - NO GO!\n" 
    and return
    unless (($params{genotype_id} or $params{genotype}) and ($params{environment_id} or $params{environment}) 
	    and ($params{phenotype_id} or $params{phenotype}) and ($params{pub_id} or $params{pub}) 
	    and $params{type_id});

  my $pc_el = $ldoc->createElement('phenstatement');

  if ($params{genotype}) {
    $params{genotype_id} = create_ch_genotype(doc => $ldoc, uniquename => $params{genotype});
    delete $params{genotype};
  }

  if ($params{environment}) {
    $params{environment_id} = create_ch_environment(doc => $ldoc, uniquename => $params{environment});
    delete $params{environment};
  }

  if ($params{phenotype}) {
    $params{phenotype_id} = create_ch_phenotype(doc => $ldoc, uniquename => $params{phenotype});
    delete $params{phenotype};
  }

  if ($params{pub}) {
    $params{pub_id} = create_ch_pub(doc => $ldoc, uniquename => $params{pub},);
    delete $params{pub};
  }

  foreach my $e (keys %params) {
    next if ($e eq 'doc');
    $pc_el->appendChild(_build_element($ldoc,$e,$params{$e}));
  }

  return $pc_el;
}

# CREATE pub or pub_id element
# params
# doc - XML::DOM::Document required
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
#        creating a new pub (i.e. not null value but not part of unique key
# type - string from pub type cv
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value
# with_id - boolean optional if 1 then pub_id element is returned
# no_lookup - boolean option if 1 then default op="lookup" attribute will not be added to element
sub create_ch_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $pubid_el = $ldoc->createElement('pub_id') if $params{with_id};

    my $pub_el = $ldoc->createElement('pub');
    $pub_el->setAttribute('id',$params{macro_id}) if $params{macro_id};  

    # add an op="lookup" attribute unless no_lookup is specified
    unless ($params{no_lookup}) {
	$pub_el->setAttribute('op','lookup');
    }

    if ($params{type}) {
      $params{type_id} = create_ch_cvterm(doc => $ldoc,
				   name => $params{type},
				   cv => 'pub type',
				  );
      delete $params{type};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'no_lookup' || $e eq 'macro_id');
	$pub_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }

    if ($pubid_el) {
	$pubid_el->appendChild($pub_el);
	return $pubid_el;
    }
    return $pub_el;
}

# CREATE pubauthor element
# params
# doc - XML::DOM::Document required
# pub_id -  macro pub id or XML::DOM pub element optional to create a freestanding pubauthor element
# pub - pub uniquename optional but required if making a freestanding element unless pub_id 
# rank - positive integer required 
# surname - string required 
# editor - boolean 't' or 'f' default = 'f' so don't pass unless you want to change
# givennames - string optional
# suffix - string optional  
sub create_ch_pubauthor {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    print "WARNING -- missing required rank (positive integer) -- NO GO!\n" 
      and return unless ($params{rank});

    print "WARNING -- no surname parameter present -- I TRUST YOU ARE NOT DOING AN INSERT\n"
      unless ($params{surname});

    ## issue some warnings if surname, givennames or suffix are too long
    print "WARNING - long surname.  Truncation to 100 characters will occur" 
      if $params{surname} and (length($params{surname}) > 100);
    print "WARNING - long givennames.  Truncation to 100 characters will occur" 
      if $params{givennames} and (length($params{givennames}) > 100);
    print "WARNING - long suffix.  Truncation to 100 characters will occur" 
      if $params{suffix} and (length($params{suffix}) > 100);

    my $pa_el = $ldoc->createElement('pubauthor');

    # deal with pub info if provided - implies a freestanding pubauthor
    if ($params{pub}) {
      $params{pub_id} = create_ch_pub(doc => $ldoc, uniquename => $params{pub});
      delete $params{pub};
    }

    foreach my $e (keys %params) {
      next if ($e eq 'doc');
      $pa_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    
    return $pa_el;
}

sub create_ch_pubprop {
    my %params = @_;
    $params{parentname} = 'pub';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE pub_dbxref element
# params
# doc - XML::DOM::Document required
# pub_id - macro pub id or XML::DOM pub element optional to create freestanding pub_dbxref
# dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
# accession - string required unless dbxref_id provided
# db_id - macro db id or XML::DOM db element required unless dbxref_id provided
# db - string name of db
# is_current - string 't' or 'f' boolean default = 't' so don't pass unless
sub create_ch_pub_dbxref {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fd_el = $ldoc->createElement('pub_dbxref');

    if ($params{pub_id}) {
	$fd_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
	delete $params{pub_id};
    }

    my $ic;
    if (exists($params{is_current})) { #assign value to a var and then remove from params
      if ($params{is_current}) {
	$ic = $params{is_current};
      } else {
	$ic = 'false';
      }
      delete $params{is_current};
    }
    # create a dbxref element if necessary
    unless ($params{dbxref_id}) {
      print "WARNING - missing required parameters, NO GO.\n" and return unless 
	($params{accession} and ($params{db_id} or $params{db}));
      if ($params{db}) {
	$params{db_id} = create_ch_db(doc => $ldoc,
				      name => $params{db},
				     );
	delete $params{db};
      }
      
	
      $params{dbxref_id} = create_ch_dbxref(%params);
    }

    $fd_el->appendChild(_build_element($ldoc,'dbxref_id',$params{dbxref_id})); #add dbxref element
    $fd_el->appendChild(_build_element($ldoc,'is_current',$ic)) if $ic;

    return $fd_el;
}


# CREATE pub_relationship element
# NOTE: now can create freestanding
# params
# doc - XML::DOM::Document required
# object_id
# subject_id
# is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
# is_subject - boolean 't'          this flag indicates if the pub info provided should be 
#                                   added in as subject or object pub
# rtype_id - macro id for a relationship type  or XML::DOM cvterm element
# rtype - string for relationship type
# NOTE: if relationship name is given will be assigned to 'pub relationship type' cv
# pub_id - macro id for a pub or XML::DOM pub element required unless uniquename provided
# uniquename - uniquename of the pub - required unless pub element provided
# type_id - macro id for a pub type or XML::DOM cvterm element for pub type
#           note: this is optional for the moment but encouraged
# type - string from pub type cv
sub create_ch_pub_relationship {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    ## pub_relationship element (will be returned)
    my $fr_el = $ldoc->createElement('pub_relationship');
    $params{relation} = $fr_el;

    # deal with creating the pub bits if present
    unless ($params{object_id} and $params{subject_id}) {
      unless ($params{pub_id}) {
	print "WARNING - no required uniquename for pub\n" and return unless $params{uniquename};
	my %pub_params = (doc => $ldoc,
			  uniquename => $params{uniquename},
			 );
	delete $params{uniquename};
	
	if ($params{type}) {
	  $params{type_id} = create_ch_cvterm(doc => $ldoc,
					      name => $params{type},
					      cv => 'pub type',
					     );
	  delete $params{type};
	}
	if ($params{type_id}) {
	  $pub_params{type_id} = $params{type_id};
	  delete $params{type_id};
	}
	$params{pub_id} = create_ch_pub(%pub_params);
      }
      $params{thingy} = $params{pub_id}; # now we have a pub element to associate as subject or object
      delete $params{pub_id};
    }
    $params{rtypecv} = 'pub relationship type' if ($params{rtype} and ! $params{rtypecv});

    return create_ch_relationship(%params);
}
*create_ch_pr = \&create_ch_pub_relationship;

# CREATE strain element
# params
# doc - XML::DOM::Document required
# uniquename - string required
# name - string optional
# organism_id - organism macro id or XML::DOM organism element required if no genus and species
# genus - string required if no organism
# species - string required if no organism
# dbxref_id - dbxref macro id or XML::DOM dbxref element required unless accession and db
# accession - string optional
# version - int optional
# db - string dbname optional
# is_obsolete - boolean 't' or 'f' default = 'f' optional
# macro_id - string optional if provide then add an ID attribute to the top level element of provided value
# with_id - boolean optional if 1 then strain_id element is returned
# no_lookup - boolean option if 1 then default op="lookup" attribute will not be added to element
sub create_ch_strain {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fid_el = $ldoc->createElement('strain_id') if $params{with_id};

    ## strain element (will be returned)
    my $f_el = $ldoc->createElement('strain');
    $f_el->setAttribute('id',$params{macro_id}) if $params{macro_id};    

    ## add an op="lookup" attribute unless no_lookup is specified
    #unless ($params{no_lookup}) {
#	$f_el->setAttribute('op','lookup');
#    }
	
    #create organism_id element if genus and species are provided
    unless ($params{organism_id}) {
	_add_orgid_param(\%params);
    }

    if ($params{accession}) {
      _add_dbxrefid_param(\%params);
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id' || $e eq 'no_lookup');
	$f_el->appendChild(_build_element($ldoc, $e,$params{$e}));
    }

    if ($fid_el) {
	$fid_el->appendChild($f_el);
	return $fid_el;
    }
    return $f_el;
}

# CREATE strain_cvterm element
# params
# doc - XML::DOM::Document required
# cvterm_id - cvterm macro id or XML::DOM cvterm element unless other cvterm bits are provided
# name - string required unless cvterm_id provided Note: a cvterm has a lookup by default cannot make a new cvterm 
#                                                        with this method
# cv_id - macro id for cv or XML::DOM cv element required if name and not cv
# cv - string for name of cv required if name and not cv_id
# pub_id - macro id for pub or XML::DOM pub element required unless pub
# pub - string = pub uniquename Note: as pub has lookup option by default can't make a new pub using this param
# is_not - optional boolean 't' or 'f' with default = 'f' so don't pass unless you know you want to change
sub create_ch_strain_cvterm {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fct_el = $ldoc->createElement('strain_cvterm');

    #create a cvterm element if necessary
    unless ($params{cvterm_id}) {
      print "WARNING -- you haven't provided params to make a cvterm, Sorry\n" and return unless ($params{name});
      my %cvtparams = (doc => $ldoc,
		       name => $params{name},
		      );
      delete $params{name};

      if ($params{cv_id}) {
	$cvtparams{cv_id} = $params{cv_id};
	delete $params{cv_id}; 
      } elsif ($params{cv}) {
	$cvtparams{cv} = $params{cv};
	delete $params{cv}; 
      } else {
	print "WARNING -- you're trying to make a cvterm without providing a cv - Sorry, NO GO\n" and return;
      }
      $params{cvterm_id} = create_ch_cvterm(%cvtparams);      
    }

    #create a pub element if necessary
    if ($params{pub}) {
      $params{pub_id} = create_ch_pub(doc => $ldoc,
				      uniquename => $params{pub},
				     );
      delete $params{pub};
    }

    foreach my $e (keys %params) {
	next if ($e eq 'doc');
	$fct_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }
    return $fct_el;
}



# CREATE strain_cvtermprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - string from  strain_cvtermprop type cv or XML::DOM cvterm element required 
#        Note: will default to making a strainprop from 'strain_cvtermprop type' cv unless 
#              cvname is provided
# cvname - string optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_strain_cvtermprop {
    my %params = @_;
    $params{parentname} = 'strain_cvterm';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}


# CREATE strain_dbxref element
# params
# doc - XML::DOM::Document required
# strain_id - macro strain id or XML::DOM strain element optionaal to create freestanding strain_dbxref
# dbxref_id - macro dbxref id or XML::DOM dbxref element - required unless accession and db provided
# accession - string required unless dbxref_id provided
# db_id - macro db id or XML::DOM db element required unless dbxref_id provided
# db - string name of db
# version - string optional
# description - string optional
# is_current - string 't' or 'f' boolean default = 't' so don't pass unless
#              this shoud be changed
sub create_ch_strain_dbxref {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fd_el = $ldoc->createElement('strain_dbxref');

    if ($params{strain_id}) {
	$fd_el->appendChild(_build_element($ldoc,'strain_id',$params{strain_id}));
	delete $params{strain_id};
    }

    my $ic;
    if (exists($params{is_current})) { #assign value to a var and then remove from params
      if ($params{is_current}) {
	$ic = $params{is_current};
      } else {
	$ic = 'false';
      }
      delete $params{is_current};
    }
    
    # create a dbxref element if necessary
    unless ($params{dbxref_id}) {
      print "WARNING - missing required parameters, NO GO.\n" and return unless 
	($params{accession} and ($params{db_id} or $params{db}));
      if ($params{db}) {
	$params{db_id} = create_ch_db(doc => $ldoc,
				      name => $params{db},
				     );
	delete $params{db};
      }
      
	
      $params{dbxref_id} = create_ch_dbxref(%params);
    }

    $fd_el->appendChild(_build_element($ldoc,'dbxref_id',$params{dbxref_id})); #add dbxref element
    $fd_el->appendChild(_build_element($ldoc,'is_current',$ic)) if $ic;

    return $fd_el;
}

# CREATE strain_feature element
# params
# doc - XML::DOM::Document required
# organism_id - macro id for organism or XML::DOM organism element
# genus - string
# species - string
# NOTE: you can use the generic paramaters in the following cases:
#       1.  you are only building either a strain or feature element and not both
#       2.  or both strain and feature have the same organism
#       otherwise use the prefixed parameters
# WARNING - if you provide both generic and prefixed parameters then the prefixed ones will be used
# strain_id - macro id for strain or XML::DOM strain element
# str_uniquename - string strain uniquename
# str_organism_id - macro id for organism or XML::DOM organism element to link to strain
# str_genus
# str_species
# str_type_id
# str_type
# feature_id - macro id for feature or XML::DOM feature element
# feat_uniquename
# feat_organism_id - macro id for organism or XML::DOM organism element to link to feature
# feat_genus
# feat_species
# feat_type_id
# feat_type
# pub_id - macroid or pub element
# pub - str pub uniquename
sub create_ch_strain_feature {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  #have to think about the checks a bit

  my $lf_el = $ldoc->createElement('strain_feature');

  # deal with feature bits
  if ($params{feat_uniquename}) {
    print "ERROR -- you don't have required parameters to make a feature, NO GO!\n" and return
      unless (($params{organism_id} or $params{feat_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{feat_genus} and $params{feat_species}))
	      and ($params{feat_type_id} or $params{feat_type}));

    unless ($params{feat_organism_id} or ($params{feat_genus} and $params{feat_species})) {
      $params{feat_organism_id} = $params{organism_id} if $params{organism_id};
      $params{feat_genus} = $params{genus} if $params{genus};
      $params{feat_species} = $params{species} if $params{species};
    }

    # gather all the feature parameters
    my %fparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /feat_(.+)/) {
	$fparams{$1} = $params{$p};
      }
    }

    $params{feature_id} = create_ch_feature(%fparams);
  }

  unless ($params{pub_id}) {
    print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
    $params{pub_id} = create_ch_pub(doc => $ldoc,
				    uniquename => $params{pub},
				   );
    delete $params{pub};
  }
  

  # likewise deal with the strain bits
  if ($params{str_uniquename}) {
    print "ERROR -- you don't have required parameters to make a strain, NO GO!\n" and return
      unless (($params{organism_id} or $params{str_organism_id} or ($params{genus} and $params{species}) 
	       or ($params{str_genus} and $params{str_species}))
	      and ($params{str_type_id} or $params{str_type}));

    unless ($params{str_organism_id} or ($params{str_genus} and $params{str_species})) {
      $params{str_organism_id} = $params{organism_id} if $params{organism_id};
      $params{str_genus} = $params{genus} if $params{genus};
      $params{str_species} = $params{species} if $params{species};
    }

    # gather all the strain parameters
    my %lparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      if ($p =~ /str_(.+)/) {
	$lparams{$1} = $params{$p};
      }
    }

    $params{strain_id} = create_ch_strain(%lparams);
  }

  # and then add the feature, strain or both to the strain_feature element
  $lf_el->appendChild(_build_element($ldoc,'strain_id',$params{strain_id})) if $params{strain_id};
  $lf_el->appendChild(_build_element($ldoc,'feature_id',$params{feature_id})) if $params{feature_id};
  $lf_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id})) if $params{pub_id};

  return $lf_el;
}


# CREATE strain_featureprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for strain_featureprop type or XML::DOM cvterm element required
# type -  string from  strain_featureprop type cv 
#        Note: will default to making a featureprop from 'strain_featureprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'strain_featureprop type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_strain_featureprop {
    my %params = @_;
    $params{parentname} = 'strain_feature';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE strain_phenotype element
# params
# doc - XML::DOM::Document required
# strain_id - optional macro id for strain or XML::DOM strain element
# phenotype_id - optional macro id for phenotype or XML::DOM phenotype element 
# uniquename - string required if no phenotype_id
# observable_id - macro id for observable or XML::DOM cvterm element optional
# attr_id - macro id for attr or XML::DOM cvterm element optional
# cvalue_id - macro id for cvalue or XML::DOM cvterm element optional
# assay_id - macro id for assay or XML::DOM cvterm element optional
# pub_id - macro id for pub or XML::DOM pub element
# pub - string uniquename for pub
# value - string optional
# macro_id - optional string to specify as ID attribute for phenotype

sub create_ch_strain_phenotype {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document

  #have to think about the checks a bit

  my $lf_el = $ldoc->createElement('strain_phenotype');

  # deal with phenotype bits
  unless ($params{phenotype_id}) {
    print "ERROR -- you don't have required parameters to make a phenotype, NO GO!\n" and return
      unless $params{uniquename};

    # gather all the feature parameters
    my %pparams = (doc => $ldoc,);
    foreach my $p (keys %params) {
      next if ($p eq 'strain_id');
      $pparams{$p} = $params{$p};
    }

    $params{phenotype_id} = create_ch_phenotype(%pparams);
  }

  unless ($params{pub_id}) {
    print "WARNING - missing required pub info, NO GO.\n" and return unless $params{pub};
    $params{pub_id} = create_ch_pub(doc => $ldoc,
				    uniquename => $params{pub},
				   );
    delete $params{pub};
  }
  

  # and then add the feature, strain or both to the strain_feature element
  $lf_el->appendChild(_build_element($ldoc,'strain_id',$params{strain_id})) if $params{strain_id};
  $lf_el->appendChild(_build_element($ldoc,'phenotype_id',$params{phenotype_id})) if $params{phenotype_id};
  $lf_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id})) if $params{pub_id};

  return $lf_el;
}

# CREATE strain_phenotypeprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id for strain_phenotypeprop type or XML::DOM cvterm element required
# type -  string from  strain_phenotypeprop type cv 
#        Note: will default to making a featureprop from 'strain_phenotypeprop type' cv unless cvname is provided
# cvname - string (probably want to pass 'strain_phenotypeprop type') optional
#          but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_strain_phenotypeprop {
    my %params = @_;
    $params{parentname} = 'strain_phenotype';
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}

# CREATE strain_pub element
# Note that this is just calling create_ch_pub setting with_id = 1 
#      and adding returned pub_id element as a child of strain_pub
#      or just appending the pub element to strain_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - macro_id for pub or XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
#        creating a new pub (i.e. not null value but not part of unique key
# type - string for type from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_strain_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('strain_pub');

    if ($params{strain_id}) {
	$fp_el->appendChild(_build_element($ldoc,'strain_id',$params{strain_id}));
	delete $params{strain_id};
    }
    
    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}

# CREATE strain_relationship element
# NOTE: this will create either a subject_id or object_id strain_relationship element
# but you have to attach this to the related strain elsewhere
# params
# doc - XML::DOM::Document required
# object _id - macro id for object strain or XML::DOM strain element
# subject_id - macro id for subject strain or XML::DOM strain element
# NOTE you can pass one or both of the above parameters with the following rules:
# if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other strain
# if both are passed then is_object, is_subject and any parameters to create a strain are ignored
# is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
# is_subject - boolean 't'          this flag indicates if the strain info provided should be 
#                                   added in as subject or object strain
# rtype_id - macro id for relationship type or XML::DOM cvterm element (Note: currently all is_relationship = '0'
# rtype - string for relationship type note: with this param  type will be assigned to relationship_type cv unless rtypecv provided
# rtypecv - string for name of cv for type of relationship - defaults to 'relationship type'
# strain_id - macro_id for a strain or XML::DOM strain element required unless minimal strain bits provided
# uniquename - string required unless strain provided
# organism_id - macro id for organism or XML::DOM organism element required unless strain or (genus & species) provided
# genus - string required unless strain or organism provided
# species - string required unless strain or organism provided
# stype_id - macro id for strain type or XML::DOM cvterm element required unless strain provided
# stype - string for strain type
sub create_ch_strain_relationship {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document
  
  ## strain_relationship element (will be returned)
  my $fr_el = $ldoc->createElement('strain_relationship');
  $params{relation} = $fr_el;
  
  # deal with creating the strain bits if present and needed
  unless ($params{object_id} and $params{subject_id}) {
    unless ($params{strain_id}) {
      unless ($params{organism_id}) {
	_add_orgid_param(\%params);
      }
      
      # before creating strain figure out which parameters
      print "WARNING - no required uniquename for strain\n" and return unless $params{uniquename};
      my %fparams = (doc => $ldoc,
		     uniquename => $params{uniquename},
		     organism_id => $params{organism_id},
		    );
      delete $params{uniquename};
      delete $params{organism_id};
      if ($params{stype_id}) { 
	$fparams{type_id} = $params{stype_id};
	delete $params{stype_id};
      } elsif ($params{stype}) {
	$fparams{type} = $params{stype};
	delete $params{stype};
      } else {
	print "WARNING -- you need to provide a strain type to make a strain!\n" and return;
      }
      $params{strain_id} = create_ch_strain(%fparams);
    } # now we have a strain element to associate as subject or object
    
    $params{thingy} = $params{strain_id};
    delete $params{strain_id};
  } 
  #$params{rtypecv} = 'relationship type' if ($params{rtype} and ! $params{rtypecv});   
  return create_ch_relationship(%params);
}
  

# CREATE strain_relationship_pub element
# # Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of strain_relationship_pub
#      or just appending the pub element to strain_relationship_pub if that is passed
# params
# doc - XML::DOM::Document required
# strain_relationship_id -optional strain_relationship XML::DOM element or macro id
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_strain_relationship_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('strain_relationship_pub');

    if ($params{strain_relationship_id}) {
	$fp_el->appendChild(_build_element($ldoc,'strain_relationship_id',$params{strain_relationship_id}));
	delete $params{strain_relationship_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}



# CREATE strain_synonym element
# params
# doc - XML::DOM::Document required
# synonym_id - macro id for synonym? or XML::DOM synonym element required unless name and type provided
# name - string required unless synonym_id element provided
# type_id - macro id for synonym type or XML::DOM cvterm element
# type - string = name from the 'synonym type' cv
# pub_id - macro id for pub or a XML::DOM pub element required
# pub - a pub uniquename (i.e. FBrf)
# synonym_sgml - string optional but if not provided then synonym_sgml = name
#               - do not provide if a synonym element is provided
# is_current - optional string = 'f' or 't' default is 't' so don't provide this param
#               unless you know you want to change the value
# is_internal - optional string = 't' default is 'f' so don't provide this param
#               unless you know you want to change the value
sub create_ch_strain_synonym {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $ls_el = $ldoc->createElement('strain_synonym');

    # this is exactly the same thing as create_ch_feature_synonym with different name
    # so call that method then change parentage?
    my $el = create_ch_feature_synonym(%params);
    my @children = $el->getChildNodes;
    $ls_el->appendChild($_) for @children;

    return $ls_el;
}


# CREATE strainprop element
# params
# doc - XML::DOM::Document required
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro id or XML::DOM cvterm element required 
#        Note: will default to making a strainprop from 'property type' cv unless cvname is provided
# type - string from property type cv
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_strainprop {
    my %params = @_;
    $params{parentname} = 'strain';
    unless ($params{type_id}) {
	$params{cvname} = 'property type' unless $params{cvname};
    }
    my $fp_el = create_ch_prop(%params);
    return $fp_el;
}


# CREATE strainprop_pub element
# Note that this is just calling create_ch_pub 
#      and adding returned pub_id element as a child of strainprop_pub
#      or just appending the pub element to strainprop_pub if that is passed
# params
# doc - XML::DOM::Document required
# pub_id - XML::DOM pub element - if this is used then pass this and doc as only params
# uniquename - string required
# type_id - macro id for pub type or XML::DOM cvterm element optional unless 
# type -  string from pub type cv 
# title - optional string
# volumetitle - optional string
# volume - optional string
# series_name - optional string
# issue - optional string
# pyear - optional string
# pages - optional string
# miniref - optional string
# is_obsolete - optional string 't' boolean value with default = 'f'
# publisher - optional string
sub create_ch_strainprop_pub {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    my $fp_el = $ldoc->createElement('strainprop_pub');

    if ($params{strainprop_id}) {
	$fp_el->appendChild(_build_element($ldoc,'strainprop_id',$params{strainprop_id}));
	delete $params{strainprop_id};
    }

    unless ($params{pub_id}) {
	$params{pub_id} = create_ch_pub(%params); #will return a pub element
    }
    $fp_el->appendChild(_build_element($ldoc,'pub_id',$params{pub_id}));
    return $fp_el;
}


# CREATE synonym or synonym_id element
# params
# doc - XML::DOM::Document required
# name - string required
# synonym_sgml - string optional but if not provided then synonym_sgml = name
# type_id - macro id or XML::DOM cvterm element
# type - string = name from the 'synonym type' cv 
# macro_id - string to assign id attribute to this element
# with_id - boolean optional if 1 then synonym_id element is returned
sub create_ch_synonym {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document

    #assign required synonym_sgml field to name if synonym_sgml not provided
    $params{synonym_sgml} = $params{name} unless ($params{synonym_sgml});
    ## synonym_id element (will be returned)
    my $synid_el = $ldoc->createElement('synonym_id') if $params{with_id};

    my $syn_el = $ldoc->createElement('synonym');
    $syn_el->setAttribute('id',$params{macro_id}) if $params{macro_id};     
    
    # set up a type_id element if necessary
    if ($params{type}) {
      my $cv_el = create_ch_cv(doc => $ldoc,
			       name => 'synonym type',
			      );
      $params{type_id} = create_ch_cvterm(doc => $ldoc,
					  name => $params{type},
					  cv_id => $cv_el,);
      delete $params{type};
    }

    foreach my $e (keys %params) {
      next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id');
      $syn_el->appendChild(_build_element($ldoc,$e,$params{$e}));
    }

    if ($synid_el) {
	$synid_el->appendChild($syn_el);
	return $synid_el;
    }
    return $syn_el;
}

# generic helper method for creating a simple table element that does not contain elements 
# that refer to other tables (i.e. no foreign key references)
# params
# doc - XML::DOM::Document required
# elname - string required - the element (table) name that you want to create
# required - array ref to list of required elements - do not need to pass doc here
# all other parameters wanted for that element
sub _create_simple_element {
  my %params = @_;
  print "ERROR -- no XML::DOM::Document specified\n" and return unless $params{doc};
  my $ldoc = $params{doc};    ## XML::DOM::Document
  print "ERROR -- you must provide the name of the elementname\n" 
    and return unless $params{elname};  
  print "ERROR -- you must provide at least one required element\n" 
    and return unless $params{required};

  my $elname = $params{elname};
  my @required = @{$params{required}};
  delete $params{elname};
  delete $params{required};

  for (@required) {
    print "ERROR: You are missing a required parameter - $_ - NO GO!\n" and return 
      unless $params{$_};
  }

  my $id_el = $ldoc->createElement("${elname}_id") if $params{with_id};


  my $el = $ldoc->createElement("$elname");
  $el->setAttribute('id',$params{macro_id}) if $params{macro_id};  
  
  foreach my $e (keys %params) {
    next if ($e eq 'doc' || $e eq 'with_id' || $e eq 'macro_id');
    $el->appendChild(_build_element($ldoc,$e,$params{$e}));
  }
  if ($id_el) {
    $id_el->appendChild($el);
    return $id_el;
  }
  return $el;
}

# generic method to deal with relationship specific paramters in one place
# doc - XML::DOM::Document required
# relation - empty *_relationship element for the appropriate table
# object _id - macro id for object or XML::DOM strain element
# subject_id - macro id for subject or XML::DOM strain element
#      NOTE you can pass one or both of the above parameters with the following rules:
#      if only one of the two are passed then the converse is_{object,subject} param is assumed for creation of other member
# if both are passed then is_object, is_subject and element parameter is ignored
# is_object - boolean 't'     Note: either is_subject OR is_object and NOT both must be passed
# is_subject - boolean 't'              this flag indicates if the element parameter is provided should be 
#                                                   as subject or object element
# thingy - an XML::DOM element or macro id for the thingy to be added as a relationship
# rtype_id - macro id for relationship type or XML::DOM cvterm element (Note: currently all is_relationship = '0'
# rtype - string for relationship type note: with this param  type will be assigned to relationship type cv unless rtypecv is provided
# rtypecv - string name of cv with relationship types
# rank - int optional default = 0 for those with this but only 2 _relationship tables (feature, strain) currently have this
# value - string optional for those with this but only 2 _relationship tables (feature, strain) currently have this
sub create_ch_relationship {
  my %params = @_;
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
  print "WARNING - the function call did not provide the top level element to build up\n" and return unless $params{relation};

  my $ldoc = $params{doc};

  my $r_el = $params{relation};

  # deal with the type
  print "WARNING - You haven't provided relationship type info -- NO GO!\n" and return unless ($params{rtype} or $params{rtype_id});
  my $rtype_el = $ldoc->createElement('type_id');
  if ($params{rtype}) {
    my $cvname = 'relationship type';
    $cvname = $params{rtypecv} if $params{rtypecv};
    my $cv_el = create_ch_cv(doc => $ldoc,
			     name => $cvname,
			    );
    my $cvterm_el = create_ch_cvterm(doc => $ldoc,
				     name => $params{rtype},
				     cv_id => $cv_el,
				     # note that we may want to add is_relationship = 1 in the future
				    );
    $rtype_el->appendChild($cvterm_el);
  } else {
    if (!ref($params{rtype_id})) {
      my $val = $ldoc->createTextNode("$params{rtype_id}");
      $rtype_el->appendChild($val);
    } else {
      $rtype_el->appendChild($params{rtype_id});
    }
  }
  $r_el->appendChild($rtype_el);
  
  # now deal with various subject/object options
  if ($params{object_id}) { 
    # create the object_id element
    $r_el->appendChild(_build_element($ldoc,'object_id',$params{object_id}));
    $params{is_subject} = 1 unless $params{subject_id};
  }
  
  if ($params{subject_id}) { 
    # create the subject_id element
    $r_el->appendChild(_build_element($ldoc,'subject_id',$params{subject_id}));
    $params{is_object} = 1 unless $params{object_id};
  }
  
  $r_el->appendChild(_build_element($ldoc,'rank',$params{rank})) if $params{rank};
  $r_el->appendChild(_build_element($ldoc,'value',$params{value})) if $params{value};

  return $r_el if (defined($params{object_id})  and defined($params{subject_id}));

  my $fr_id;
  if ($params{is_object}) {
    $fr_id = $ldoc->createElement('object_id');
  } else {
    $fr_id = $ldoc->createElement('subject_id');
  }

  if (!ref($params{thingy})) {
    my $val = $ldoc->createTextNode("$params{thingy}");
    $fr_id->appendChild($val);
  } else {
    $fr_id->appendChild($params{thingy});
  }
  
  $r_el->appendChild($fr_id);

 
  return $r_el;
}

# generic method for creating a prop element
# params
# doc - XML::DOM::Document required
# parentname - string required - the parent element name that you want to add prop to eg. feature
#              to make a featureprop or feature_cvterm to make feature_cvtermprop
# parent_id - optional macro_id or XML::DOM element for an table_id for the prop to make a standalone prop
#             NOTE that this param name should be eg. feature_id for featureprop or pub_id for pubprop etc.
# value - string - not strictly required but if you don't provide this then not much point
# type_id - macro_id for a property type or XML::DOM cvterm element required
# type - string from a property type cv
#        Note: will default to making a type of from 'tablenameprop type' cv unless cvname is provided
# WARNING: as property type cv names are not consistent SAFEST to provide cvname
# cvname - string 'optional' but see above for type and do not provide if passing a cvterm element
# rank - integer with a default of zero so don't use unless you want a rank other than 0
sub create_ch_prop {
    my %params = @_;
    print "WARNING -- no XML::DOM::Document specified\n" and return unless $params{doc};
    my $ldoc = $params{doc};    ## XML::DOM::Document
    print "WARNING -- you must provide the name of the element to which this prop will be added\n" 
      and return unless $params{parentname};
    my $tbl = $params{parentname};
    print "WARNING -- you haven't provided the required information for type -- NO GO!\n" 
      and return unless ($params{type} or $params{type_id});

    my $p_el = $ldoc->createElement("${tbl}prop");

    # check to see if making a stand alone prop element
    if ($params{"${tbl}_id"}) {
      $p_el->appendChild(_build_element($ldoc,"${tbl}_id",$params{"${tbl}_id"}));
      delete $params{"${tbl}_id"};
    }

    foreach my $e (keys %params) {
	next unless ($e eq 'value' || $e eq 'type_id' || $e eq 'type');
	if ($e eq 'value') {
	    #my $val = $params{$e};
	    ## Strip non-ascii characters from the featureprop value
	    #$val =~ s/[\x80-\xff]//g;
	    ## Escape single-quotes
	    #$val =~ s/\'/\\\'/g;
	    #$params{$e} = $val;
	  unless (!defined($params{value}) or $params{value} eq '') {
	    $p_el->appendChild(_build_element($ldoc,$e,$params{$e}));
	  }
	} elsif ($e eq 'type_id') {
	  $p_el->appendChild(_build_element($ldoc,$e,$params{$e}));
	} else {
	  my $cv = "${tbl}prop type";
	  $cv = $params{cvname} if $params{cvname};
	  my $cv_el =  create_ch_cv(doc => $ldoc,
				    name => $cv,
				   );
	  my $ct_el = (create_ch_cvterm(doc => $ldoc,
					name => $params{$e},
					cv_id => $cv_el,));
	  $p_el->appendChild(_build_element($ldoc,'type_id',$ct_el));
	}
    }

    #add the required rank element 
    my $rank = '0';
    $rank = $params{rank} if $params{rank};
    $p_el->appendChild(_build_element($ldoc,'rank',$rank));

    return $p_el;
}

# helper method to deal with common dbxref processing
sub _add_dbxrefid_param {
  my $params = shift;
  my $doc = $params->{doc};
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $doc;
  print "WARNING -- no dbxref info provided -- NO GO!\n" and return
    unless ($params->{accession} and ($params->{db_id} or $params->{db}));
  my %dxp;
  unless ($params->{db_id}) {
    my $db_el = create_ch_db(doc => $doc, name => $params->{db});
    $dxp{db_id} = $db_el;
    delete $params->{db};
  }
  $dxp{accession} = $params->{accession};
  delete $params->{accession};
  if ($params->{version}) {
    $dxp{version} = $params->{version};
    delete $params->{version};
  }
  if ($params->{description}) {
    $dxp{description} = $params->{description};
    delete $params->{description};
  }
  if ($params->{url}) {
    $dxp{url} = $params->{url};
    delete $params->{url};
  }
  $dxp{doc} = $doc;
  $params->{dbxref_id} = create_ch_dbxref(%dxp);
  delete $params->{db_id} if $params->{db_id};
  undef %dxp;
  return;
}

# helper method to deal with common organism processing
sub _add_orgid_param {
  my $params = shift;
  my $doc = $params->{doc};
  print "WARNING -- no XML::DOM::Document specified\n" and return unless $doc;
  print "WARNING -- no organism info provided -- NO GO!\n" and return
	    unless ($params->{genus} and $params->{species});
  $params->{organism_id} = create_ch_organism(doc => $doc,
						  genus => $params->{genus},
						  species => $params->{species},
					     ) ;
  delete $params->{genus}; delete $params->{species};
  return;
}

# helper method to build up elements
sub _build_element {
    my $doc = shift;
    my $ename = shift;
    my $eval = shift;

    my $el = $doc->createElement("$ename");
    my $val;
    if (!ref($eval)) {
	$val = $doc->createTextNode("$eval");
    } else {
	$val = $eval;
    }
    $el->appendChild($val);
    return $el;
}



1;
