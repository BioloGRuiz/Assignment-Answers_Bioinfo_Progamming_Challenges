### Assignment 4_Bioinformatics Programming Challenges ###

## Done by Sergio Gonzalez Ruiz

require 'bio'
#require 'ruby-progressbar'

arabidopsis_th = ARGV[0]
s_pombe = ARGV[1]
unless arabidopsis_th && s_pombe 
  abort "ERROR: try running this using the following command:\n ruby assignment4_test.rb ./BLAST_databases/arabidopsis_th.fa ./BLAST_databases/s_pombe.fa"
end

# Bio::FlatFile allows us to read a biological data file. It can automatically detect data format, and users do not need to tell the class what the data is.
# documentation --> http://bioruby.org/rdoc/Bio/FlatFile.html
arabth_flat = Bio::FlatFile.auto(arabidopsis_th)
s_pombe_flat = Bio::FlatFile.auto(s_pombe)

# Now, we will obtain the type of the sequences given: Amino Acid or Nucleic Acid, and return a new Bio::Sequence object wrapping the sequences of the guessed type (either Bio::Sequence::AA or Bio::Sequence::NA)
# documentation --> http://bioruby.org/rdoc/Bio/Sequence.html#method-c-auto
def database_type(flatfasta)
 input = Bio::Sequence.auto(flatfasta.next_entry.to_s).guess
  if input == Bio::Sequence::NA
    type = 'nucl'
  elsif input == Bio::Sequence::AA
    type = 'prot'
  end
  return type
end

# Now, we will create a function which will obtain the databases that are needed to do the reciprocal Blasts:
def make_blast_db(flatfasta,fasta)
 type = database_type(flatfasta)
 system("makeblastdb -in #{fasta} -dbtype '#{type}' -out $(dirname #{fasta})/$(basename #{fasta} .fa)") # This command needs an input fasta file and the type of sequences in the file (provided above)
end

# The following function will perform the BLAST analysis --> It performs the best reciprocal hit analysis between the two databases previously created
def blast(db,query,type)
 evalue = '-e 1e-6' # We will define an e-value threshold of 1×10−6 to find Reciprocal Best Hits (RBH), as indicated in the bibliography:
 # e-value threshold references:
 ## https://doi.org/10.1093/bioinformatics/btm585
 ## https://dx.doi.org/10.1371%2Fjournal.pone.0101850
 ## https://doi.org/10.1186/s12864-020-07132-6
 # Now we will create the factory (which will execute homology search)
 factory = Bio::Blast.local("#{type}","#{File.dirname(db)}/#{File.basename(db,".fa")}","-F 'm S' #{evalue}")
 report = factory.query(query)
 # Searching the query with the definition of the first hit
 if report.hits[0]
  return report.hits[0].definition.split("|")[0].strip
 end
end

# Now, we will simply add the types needed for this assignment:
def blast_type(flat,queries)
  if database_type(flat) == 'prot' and database_type(queries) == 'nucl'
    type='blastx'
  elsif database_type(flat) == 'nucl' and database_type(queries) == 'prot'
    type='tblastn'
  end
 return type
end


def get_best_reciprocal_hits(db,flat,queries)
 best_reciprocal_hits = Hash.new
 type = blast_type(flat,queries)
 count = 0  # UNCOMMENT THIS FOR TESTING
  queries.each_entry do |query|
   puts "Using query: #{query.entry_id}" ## Just to see that the script is working ## USE THIS JUST FOR TESTING WITH FEW ITERATIONS!!
   best_reciprocal_hits[query.entry_id] = blast(db,query,type)
   break if count == 100 # Just to make it shorter, it will stop when reaching 100 queries. UNCOMMENT THIS FOR TESTING
   count += 1  # UNCOMMENT THIS FOR TESTING
  end
 return best_reciprocal_hits
end

# With this function we will search for reciprocal best hits between the 2 species, having into account the BLAST analysis previously performed
def get_orthologues(db1,db2,flat1,flat2)
 orthologues = Hash.new
 puts "Searching for reciprocal best hits (this will take a while...)"
 best_hits = get_best_reciprocal_hits(db1,flat1,flat2)
 flat1.each_entry do |input|
  next unless best_hits.value?(input.entry_id)
  reciprocal_hits=blast(db2,input,blast_type(flat2,flat1))
   if best_hits[reciprocal_hits] == input.entry_id
    orthologues[reciprocal_hits] = input.entry_id
   end
  end
  return orthologues 
end

# Finally, calling the different functions:
make_blast_db(arabth_flat,arabidopsis_th)
make_blast_db(s_pombe_flat,s_pombe)
get_best_reciprocal_hits(arabidopsis_th,arabth_flat,s_pombe_flat) # JUST FOR TESTING
orthologues = get_orthologues(arabidopsis_th,s_pombe,arabth_flat,s_pombe_flat)

# Creating the .txt file, in which we will introduce the information obtained from the putative orthologues found:
count = 0
File.open('orthologues_test.txt', 'w+') do |orth|
 orth.puts "Possible orthologues combinations found (from S. pombe and A. thaliana, respectively):\n\n"
 orthologues.each do |pombe,thaliana|
  orth.puts "\tPair #{count + 1}: #{pombe} and #{thaliana}"
  break if count == 100 # Just to make it shorter, it will stop when reaching 100 queries. # UNCOMMENT THIS FOR TESTING
  count+=1
 end
 orth.puts "\nThe total number of putative orthologues detected is #{count}"
end
