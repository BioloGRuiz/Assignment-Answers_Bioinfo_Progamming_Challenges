require 'rest-client'
require 'json'
require './utils.rb'
require  './gene.rb'
require  './interactionNetwork.rb'

def prueba_fetch()
  #res = Utils.fetch('http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=At3g54340')
  res = Utils.fetch("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/At3g54340?species:3702?format=xml25")

  puts res
  
  if res  # res is either the response object (RestClient::Response), or false, so you can test it with 'if'
    body = res.body  # get the "body" of the response
    #headers = res.headers  # get other details about the HTTP message itself
    #puts body
    #puts res
  else
    puts "the Web call failed - see STDERR for details..."
  end
end


def read_txt(path)
  #How to read the txt: https://www.rubyguides.com/2015/05/working-with-files-ruby/
  genes_list = File.read(path).split
  #How to put caps the whole array: https://stackoverflow.com/questions/11402362/how-can-i-uppercase-each-element-of-an-array/11402608
  genes_list.map!(&:upcase)
  return genes_list
end

def prueba_gene()  
  genes_list_names = read_txt("ArabidopsisSubNetwork_GeneList.txt")
  genes_list = Array.new
  genes_list_names.each do |gene_id|
      gene = Gene.new(id: gene_id)
      genes_list << gene
  end
  
  genes_list.each do |gene|
    puts gene.gene_id
  end
end

#################
################# HASTA AQUI LAS PRUEBAS
#################


depth = 3

#cargar genes
genes_list_names = read_txt("Corto_Pruebas.txt")

genes_list = Array.new
genes_list_names.each do |gene_id|
    gene = Gene.new(id: gene_id, padre: nil, depth: depth)
    genes_list << gene
end

net_list = Array.new
genes_list.each do |gene|
  new_net = InteractionNetwork.new(gene_id: gene.gene_id, depth: depth)
  net_list << new_net
  
  puts new_net.network
  puts ''
end



