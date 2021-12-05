## Assignment done by Sergio Gonzalez Ruiz
require 'bio'             # Calling the "BioRuby library"
require 'rest-client'     # Calling the HTTP and REST client for Ruby

# We are going to call the list of 167 genes --> "ArabidopsisSubNetwork_GeneList.txt" using ARGV[0]

#---FETCH----#  Obtained from the slides of this subject (why do we need to reinvent the wheel? :P)
def self.fetch(url, headers = {accept: "*/*"}, user = "", pass="")
  response = RestClient::Request.execute({
    method: :get,
    url: url.to_s,
    user: user,
    password: pass,
    headers: headers})
  return response
    
  rescue RestClient::ExceptionWithResponse => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue RestClient::Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
end

## 1) We will retrieve the sequences of the Arabidopsis genes, using BioRuby:

def get_embl(file)
  loci_codes=Hash.new # Saving the AGI Loci codes and their EMBL entries in an instance variable hash
    File.open(file).each do |gene|
      gene.strip! # Removing leading and trailing whitespace from str. Returns nil if str was not altered --> https://apidock.com/ruby/String/strip%21
      response = fetch("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene}&style=raw") # Database selected --> Ensembl Genomes Gene (ensemblgenomesgene)
      if response
        embl=Bio::EMBL.new(response.body) # Here we will apply the BioRuby Class "Bio::EMBL", which enable us to represent the database records selected
        loci_codes[gene]=embl # Saving the AGI Loci codes and their EMBL entries in the hash previously defined
      end
    end
  return loci_codes
end

# The Bio::EMBL object finds every 'CTTCTT' sequence in the exons of each gene.
# We will add the features selected to a new Bio::Sequence object, which will contain all the previous information of the Bio::EMBL object


def add_features(pos_list, strand, bioseq) # Adding new features to the Sequence entries
  #pos_list --> list in which we introduced the coordinates of every CTTCTT sequence
  #strand --> DNA strand orientation. Positive --> '+', negative --> '-' (complementary). BioRuby reads and presents sequences in 5' → 3'
  #bioseq --> Bio::Sequence object, made from a Bio::EMBL object to be able to add features to it
  pos_list.each do |pos|
    ft=Bio::Feature.new('CTTCTT_repetition',pos) # just creating a Bio::Feature object --> http://bioruby.org/rdoc/Bio/Feature.html
    ft.append(Bio::Feature::Qualifier.new('repeat_motif','cttctt')) # here we add some qualifier-value pairs for sequence features --> http://bioruby.org/rdoc/Bio/Feature/Qualifier.html
    ft.append(Bio::Feature::Qualifier.new('function','insertion site')) # here we add some qualifier-value pairs for sequence features --> http://bioruby.org/rdoc/Bio/Feature/Qualifier.html
    if strand == 1
      ft.append(Bio::Feature::Qualifier.new('strand', '+')) # here we add some qualifier-value pairs for sequence features
    elsif strand == -1
      ft.append(Bio::Feature::Qualifier.new('strand', '-')) # here we add some qualifier-value pairs for sequence features
    end
    bioseq.features << ft # just appending all features obtained to the Bio::Sequence object
  end
end

def load_new_data(gene_list)
  new_gene_list = Hash.new # creating a Hash in which we will store each of the sequences found in the exons, along with their features
  forward_seq=(Bio::Sequence::NA.new("CTTCTT")).to_re # Regular expression used for the + strand
  reverse_seq=(Bio::Sequence::NA.new("AAGAAG")).to_re # Regular expression used for the - strand --> Complementary
  
  gene_list.each do |code, embl|
  
    bioseq = embl.to_biosequence # just turning the Bio::EMBL into Bio::Sequence so that we can add features
    positions_added_f = [] # In this local variable we will store the different added positions of each gene, for the + strand
    positions_added_r = [] # In this local variable we will store the different added positions of each gene, for the - strand --> Complementary
    embl.features do |feature|
      
      if feature.feature == "exon" # retrieving the location and NA sequence of every exon found
        feature.locations.each do |location|  
          exon_seq = embl.seq[location.from..location.to]
          
          if exon_seq == nil  # If it is empty, go next
            next
          end
          if location.strand == 1   # If the found sequence is located in the primary string (+ strand)
            if exon_seq.match(forward_seq) # if any exon contains the NA sequence "CTTCTT":
              # retrieve the positions for each of the "CTTCTT" sequences found.
              # matching the regex and getting the position in the string of the match --> https://coderedirect.com/questions/175103/ruby-regex-match-and-get-positions-of
              # as the sequences can be OVERLAPPED (for example, "CTTCTTCTT"), we need to use something like a lookahead assertion --> https://stackoverflow.com/questions/11430863/how-to-find-overlapping-matches-with-a-regexp
              # with the Ruby Map Method (".map"), we add to each of the positions of these found sequences the respective positions of their respective exons --> https://www.rubyguides.com/2018/10/ruby-map-method/
              positions_f = exon_seq.enum_for(:scan, /(?=(cttctt))/i).map { [location.from + Regexp.last_match.begin(0),location.from +  Regexp.last_match.begin(0) + 5 ]}
              
              #positions_f.each do |pos|  --> just a loop for printing in the terminal the sequences found in the + strand
              #  puts(exon_seq[pos[0]..pos[1]])
              #end
              
              positions_f = positions_f.map {|pos| pos.join('..')} # with the Ruby Map Method (".map"), we join the start and end points of each of the positions of these found sequences --> https://www.rubyguides.com/2018/10/ruby-map-method/
              positions_added_f |= positions_f # adding positions without duplicates --> https://stackoverflow.com/questions/8569039/ruby-assignment-operator
            end
          end
          if location.strand == -1  # If the found sequence is located in the complementary string (- strand)
            if exon_seq.match(reverse_seq) # if any exon contains the NA sequence "AAGAAG":
              # retrieve the positions for each of the "AAGAAG" sequences found.
              # matching the regex and getting the position in the string of the match --> https://coderedirect.com/questions/175103/ruby-regex-match-and-get-positions-of
              # as the sequences can be OVERLAPPED (for example, "AAGAAGAAG"), we need to use something like a lookahead assertion --> https://stackoverflow.com/questions/11430863/how-to-find-overlapping-matches-with-a-regexp
              # with the Ruby Map Method (".map"), we add to each of the positions of these found sequences the respective positions of their respective exons --> https://www.rubyguides.com/2018/10/ruby-map-method/
              positions_r = exon_seq.enum_for(:scan, /(?=(aagaag))/i).map { [location.from +  Regexp.last_match.begin(0), location.from +  Regexp.last_match.begin(0) + 5 ]}

              #positions_r.each do |pos| --> just a loop for printing in the terminal the sequences found in the - strand
              #  puts(exon_seq[pos[0]..pos[1]])
              #end
              positions_r = positions_r.map {|pos| pos.join('..')} # with the Ruby Map Method (".map"), we join the start and end points of each of the positions of these found sequences --> https://www.rubyguides.com/2018/10/ruby-map-method/
              positions_added_r |= positions_r # adding positions without duplicates --> https://stackoverflow.com/questions/8569039/ruby-assignment-operator
            end
          end
        end
      end
    end
    # Now we will add the different positions into the object --> bioseq
    add_features(positions_added_f, 1, bioseq) # Positions for the forward sequences
    add_features(positions_added_r, -1, bioseq) # Positions for the reverse sequences
    # IF BOTH ARE NOT EMPTY, ADD THEM TO THE NEW LIST
    if !(positions_added_r.empty? && positions_added_f.empty?) 
        new_gene_list[code] = bioseq unless new_gene_list.keys.include?(code) # UNLESS THEY ARE ALREADY FROM THE SAME CODE (gene)
    end
  end
  return new_gene_list
end


def report_genes_gff3(new_gene_list)  # Whith this function we will create a GFF3-formatted file of the different CTTCTT features obtained, along with the coordinates of every CTTCTT sequence obtained
  source="BioRuby"
  type="direct_repeat"
  score="."
  phase="."
  File.open('CTTCTT_genes_report_1.gff3', 'w+') do |gen| # Just writing the file in .gff3 format
    gen.puts("##gff-version 3") # The first line of a GFF3 file must be a comment that identifies the version --> http://www.ensembl.org/info/website/upload/gff3.html
    #gen.puts("GFF3 REPORT: FEATURES OBTAINED FOR THE SEQUENCE 'CTTCTT'")
    #gen.puts("AGI locus code\tSource\tType\tStartPos\tEndpos\tScore\tStrand\tPhase\tAttributes") # printing headers
    new_gene_list.each do |code, bioseq|
      contador = 0
      bioseq.features.each do |feature|
        if feature.feature == 'CTTCTT_repetition'
            contador+=1
            pos=feature.locations.first # getting the first location object
            strand=feature.assoc['strand'] # geting strand qualifiers
            attributes="ID=CTTCTT_insertional_repeat_#{code.upcase}_#{contador};" # adding different attributes for each features
            gen.puts("#{code.upcase}\t#{source}\t#{type}\t#{pos.from}\t#{pos.to}\t#{score}\t#{strand}\t#{phase}\t#{attributes}") # printing all different features obtained --> http://bioruby.org/rdoc/Bio/GFF.html
        end
      end
    end
  end
end

def report_chromosome_gff3(new_gene_list) # Whith this function we will create a GFF3-formatted file of the different CTTCTT features obtained, along with the full chromosome coordinates in which the CTTCTT regions appear
  source="BioRuby"
  type="direct_repeat"
  score="."
  phase="."
  File.open('CTTCTT_chromosomes_report_1.gff3', 'w+') do |chr| # Just writing the file in .gff3 format
    chr.puts("##gff-version 3") # The first line of a GFF3 file must be a comment that identifies the version --> http://www.ensembl.org/info/website/upload/gff3.html
    #chr.puts("GFF3 REPORT: FEATURES OBTAINED FOR THE SEQUENCE 'CTTCTT'")
    #chr.puts("Chromosome seqid\tSource\tType\tStartPos\tEndpos\tScore\tStrand\tPhase\tAttributes") # printing headers
    new_gene_list.each do |code, bioseq|
      chr_coords=bioseq.primary_accession.split(":")[3] # HERE WE SELECT THE BEGINING POSITION OF THE CHROMOSOME
      seqid="Chr#{bioseq.primary_accession.split(":")[2]}" # SELECTING THE CHROMOSOME NUMBER
      contador = 0
      bioseq.features.each do |feature|
        if feature.feature == 'CTTCTT_repetition'
            contador+=1
            pos=feature.locations.first # getting the first location object
            strand=feature.assoc['strand'] # geting strand qualifiers
            attributes="ID=CTTCTT_insertional_repeat_#{code.upcase}_#{contador};" # adding different attributes for each features
            
            first=chr_coords.to_i+pos.from # creating the locations relative to the chromosome beginning position
            last=chr_coords.to_i+pos.to    # creating the locations relative to the chromosome beginning position
            
            chr.puts("#{seqid}\t#{source}\t#{type}\t#{first}\t#{last}\t#{score}\t#{strand}\t#{phase}\t#{attributes}") # printing all different features obtained (http://www.ensembl.org/info/website/upload/gff3.html)
        end
      end
    end
  end
end

def noreps_report(gene_list,new_gene_list) # now we will create another report, in which the loci with no CTTCTT repeats are listed
  contador=0
  File.open('CTTCTT_noRepeats_report_1.txt', 'w+') do |norep|
    norep.puts("Here is presented the list of loci which contain no CTTCTT repeats:")
    gene_list.each do |a,b|
      unless new_gene_list.keys.include?(a) # the loci codes that are not in the new_gene_list variable are the ones for which the repeat has not been found
        contador+=1
        norep.puts("\t#{contador} : #{a.upcase}")
      end
    end
  end
end

#####################
#####################
#####################

# Final results
puts "Here we go, relax and go for a coffee!"
gene_list = get_embl(ARGV[0])

puts "Searching for repeats and adding features..."
new_gene_list = load_new_data(gene_list)

puts "Writing a GFF3 file with gene coordinates and features..."
report_genes_gff3(new_gene_list)

puts "Writing a GFF3 file with chromosomes coordinates and features..."
report_chromosome_gff3(new_gene_list)

puts "Writing a report with the IDs of all the genes in the initial gene list that don't contain 'CTTCTT' in any of their exons..."
noreps_report(gene_list,new_gene_list)

