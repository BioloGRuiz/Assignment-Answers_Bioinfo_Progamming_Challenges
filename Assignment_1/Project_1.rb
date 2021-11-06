class Seed_stock

  attr_accessor :seed_stock_number
  attr_accessor :mutant_gene_id
  attr_accessor :last_planted
  attr_accessor :storage
  attr_accessor :grams_remaining
  attr_accessor :gene
  
  
  def initialize (params = {}) # get a name from the "new" call, or set a default
    @seed_stock_number = params.fetch(:seed_stock_number, "00000") 
    @mutant_gene_id = params.fetch(:mutant_gene_id, "AT0G00000")
    @last_planted = params.fetch(:last_planted, "00/00/0000")
    @storage = params.fetch(:storage, [])
    @grams_remaining = params.fetch(:grams_remaining, nil).to_i # changing the original type from string to integer
    @gene = params.fetch(:gene, nil)
    
  end

  
  def planting_seeds (n)
    @grams_remaining -= n
    if @grams_remaining <= 0
      puts "WARNING: we have run out of Seed Stock #{@seed_stock_number}."
      @grams_remaining = 0
    end
    @last_planted = DateTime.now.strftime "%d/%m/%Y" # Set today as the last planted date
  end
  
end

class Gene_information

  attr_accessor :gene_id
  attr_accessor :gene_name
  attr_accessor :mutant_phenotype
  attr_accessor :linked_genes
  
  
  def initialize (params = {}) # get a name from the "new" call, or set a default
    @gene_id = params.fetch(:gene_id, "AT0G00000") 
    @gene_name = params.fetch(:gene_name, "any gene")
    @mutant_phenotype = params.fetch(:mutant_phenotype, [])
    @linked_genes = Hash.new # here we will introduce the genes that are linked
  end
  
  def add_linked_gene(gene, chisq)
    @linked_genes[gene] = chisq  
  end
end

class Hybrid_cross

  attr_accessor :p1 #parent gene object
  attr_accessor :p2 #parent gene object
  attr_accessor :f2_wild #offspring
  attr_accessor :f2_p1 #offspring
  attr_accessor :f2_p2 #offspring
  attr_accessor :f2_p1p2 #offspring
  attr_accessor :chi_sq 
  
  def initialize (params = {}) # get a name from the "new" call, or set a default
    @p1 = params.fetch(:p1, "00000") 
    @p2 = params.fetch(:p2, "00000")
    @f2_wild = params.fetch(:f2_wild, "000").to_f  ##The values are converted to decimals in order to calculate the chi square
    @f2_p1 = params.fetch(:f2_p1, "00").to_f
    @f2_p2 = params.fetch(:f2_p2, "00").to_f
    @f2_p1p2 = params.fetch(:f2_p1p2, "00").to_f
    
    # Calculating chi-square
    summation = @f2_wild + @f2_p1 + @f2_p2 + @f2_p1p2 
    @chi_sq = ((@f2_wild - summation * 9/16)**2) / (summation * 9/16) +
              ((@f2_p1 - summation * 3/16)**2) / (summation * 3/16) +
              ((@f2_p2 - summation * 3/16)**2) / (summation * 3/16) +
              ((@f2_p1p2 - summation * 1/16)**2) /(summation * 1/16)
    if @chi_sq > 7.815 # chi square value for a p value < 0.05 with 3 degrees of freedom (4 offspring --> n-1). Source: https://www.mun.ca/biology/scarr/4250_Chi-square_critical_values.html
      puts "Recording: #{@p1.gene.gene_name} is genetically linked to #{@p2.gene.gene_name} with chisquare score #{@chi_sq}"
      @p1.gene.add_linked_gene(@p2.gene, @chi_sq)
      @p2.gene.add_linked_gene(@p1.gene, @chi_sq) 
      
    end
    
  end
  
end

class Database
  attr_accessor :gene_hash
  attr_accessor :stock_hash
  attr_accessor :cross_list
  
  def initialize(params = {})
    #primero, guardar los datos
    @gene_hash = Hash.new
    @stock_hash = Hash.new
    @cross_list = []
    name_file_gen = params.fetch(:genefile)
    name_file_stock = params.fetch(:stockfile)
    name_file_cross = params.fetch(:crossfile)

    self.load_gene_data(name_file_gen)
    self.load_stock_data(name_file_stock)
    self.load_cross_data(name_file_cross)
  end

  
  def load_gene_data(file_name)
    gene_table = CSV.read(file_name, { :col_sep => "\t"})
    header_gene = gene_table.shift # link to .shift: https://stackoverflow.com/questions/11740439/how-can-i-skip-the-header-row-when-reading-a-csv-in-ruby/11740635
    gene_table.each do |gene|
      tem_gene = Gene_information.new(
        gene_id: gene[0],
        gene_name: gene[1],
        mutant_phenotype: gene[2])
      @gene_hash[gene[0]] = tem_gene
    end
  end

  
  def load_stock_data(file_name)
    stock_table = CSV.read(file_name, { :col_sep => "\t"})
    header_stock = stock_table.shift # link to .shift: https://stackoverflow.com/questions/11740439/how-can-i-skip-the-header-row-when-reading-a-csv-in-ruby/11740635
    stock_table.each do |stock|
      tem_stock = Seed_stock.new(
        seed_stock_number: stock[0],
        mutant_gene_id: stock[1],
        last_planted: stock[2],
        storage:stock[3],
        grams_remaining: stock[4],
        gene: gene_hash[stock[1]])
      @stock_hash[stock[0]] = tem_stock
    end
  end

  def load_cross_data(file_name)
    cross_table = CSV.read(file_name, { :col_sep => "\t"})
    header_cross = cross_table.shift # link to .shift: https://stackoverflow.com/questions/11740439/how-can-i-skip-the-header-row-when-reading-a-csv-in-ruby/11740635
    cross_table.each do |cross|
      @cross_list << Hybrid_cross.new(
        p1: stock_hash[cross[0]],
        p2: stock_hash[cross[1]],
        f2_wild: cross[2],
        f2_p1: cross[3],
        f2_p2: cross[4],
        f2_p1p2: cross[5])
      
    end
    
  end
  
end

######################################################################
######################################################################
######################################################################

require 'csv'


fichero_gene, fichero_stock, fichero_cross, newfile = ARGV


database = Database.new({
  :genefile => fichero_gene,
  :stockfile => fichero_stock,
  :crossfile => fichero_cross,
})

puts()
database.stock_hash.each do |stock|
  stock[1].planting_seeds(7)
end

CSV.open(newfile, "w") do |csv| ##source --> https://medium.com/@ali_schlereth/working-with-csvs-in-ruby-43005e566901
  database.stock_hash.each do |stock|
    csv << 'hola'
  end
  
end
