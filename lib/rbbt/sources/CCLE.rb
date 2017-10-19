require 'rbbt-util'
require 'rbbt/resource'

module CCLE
  extend Resource
  self.subdir = 'share/databases/CCLE'

  def self.organism(org="Hsa")
    Organism.default_code(org)
  end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  CCLE.claim CCLE.gene_CNV, :proc do  |filename|
    url = "https://data.broadinstitute.org/ccle_legacy_data/dna_copy_number/CCLE_copynumber_byGene_2013-12-03.txt.gz"
    tsv = TSV.open(url, :wget_options => {"--no-check-certificate" => true}, :gzip => true, :type => :list, :header_hash => '')
    fields = tsv.fields
    tsv.key_field = "Entrez Gene ID"
    good_fields = fields[4..-1]
    
    tsv.slice(good_fields).to_s
  end

  CCLE.claim CCLE.gene_expression, :proc do  |filename|
    url = "https://data.broadinstitute.org/ccle_legacy_data/mRNA_expression/CCLE_Expression_Entrez_2012-09-29.gct"
    io = Open.open(url, :wget_options => {"--no-check-certificate" => true})
    io2 = CMD.cmd('tail -n +3 | cut -f 2-', :in => io, :pipe => true)
    tsv = TSV.open(io2, :type => :list, :merge => true, :header_hash => '', :cast => :to_f)
    tsv.key_field = "Associated Gene Name"
    tsv.to_s
  end

  CCLE.claim CCLE.cell_lines, :proc do  |filename|
    url = "https://data.broadinstitute.org/ccle_legacy_data/cell_line_annotations/CCLE_sample_info_file_2012-10-18.txt"
    io = Open.open(url, :wget_options => {"--no-check-certificate" => true})
    tsv = TSV.open(io, :type => :list, :header_hash => '')
    tsv.to_s
  end

  CCLE.claim CCLE.oncomap_maf, :proc do  |filename|
    url = "https://data.broadinstitute.org/ccle_legacy_data/oncomap_mutations/CCLE_Oncomap3_2012-04-09.maf"
    io = Open.open(url, :wget_options => {"--no-check-certificate" => true})
    tsv = TSV.open(io, :type => :double, :header_hash => '', :merge => true)
    tsv.key_field = "Associated Gene Name"
    tsv = tsv.reorder "Tumor_Sample_Barcode"
    tsv.add_field "Genomic Mutation" do |sample, values|
      mutations = Misc.zip_fields(values.values_at("Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")).collect do |chr,start,ref,alt1,alt2|
        pos, muts = Misc.correct_mutation(start, ref, [alt1,alt2].compact.uniq * ",")
        muts.collect{|mut| [chr,pos,mut] * ":" }.flatten.uniq
      end
      mutations
    end
    tsv.to_s
  end

  CCLE.claim CCLE.drugs, :proc do  |filename|
    url = "https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_profiling_2012.02.20.csv"
    io = Open.open(url, :wget_options => {"--no-check-certificate" => true})
    tsv = TSV.open(io, :type => :list, :sep => ',', :header_hash => '', :fix => Proc.new{|l| l.gsub(/"[^"]+"/){|p| p.gsub(', ',"|").gsub('"','')}} )
    tsv.to_s
  end

  CCLE.claim CCLE.drug_profiles, :proc do  |filename|
    url = "https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv"

    io = Open.open(url, :wget_options => {"--no-check-certificate" => true})
    tsv = TSV.open(io, :type => :double, :merge => true, :sep => ',', :header_hash => '', :fix => Proc.new{|l| l.gsub(/"[^"]+"/){|p| p.gsub(',',"|").gsub('"','')}} )
    tsv.to_s
  end

  CCLE.claim CCLE.cell_line_binary, :proc do  |filename|
    url = "https://data.broadinstitute.org/ccle_legacy_data/binary_calls_for_copy_number_and_mutation_data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct"
    io = Open.open(url, :wget_options => {"--no-check-certificate" => true})
    io2 = CMD.cmd('tail -n +3 | cut -f 2-', :in => io, :pipe => true)
    tsv = TSV.open(io2, :type => :list, :header_hash => '', :cast => :to_i)
    tsv.to_s
  end

end

Log.with_severity 0 do

  Log.tsv CCLE.drug_profiles.produce(true).tsv 

end if __FILE__ == $0

