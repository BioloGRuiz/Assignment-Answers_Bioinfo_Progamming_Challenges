require './utils.rb'
require  './gene.rb'
require  './interactionNetwork.rb'

#gene_list_txt, final_report = ARGV
#gene_list_txt, final_report =
 # unless gene_list_txt && report 
  #  abort "run this using the command\n main.rb ArabidopsisSubNetwork_GeneList.txt final_report.txt"
  #end

def write_report(final_report, all_networks, gene_txt)
  

  # REPORT CREATION
  
  File.open(final_report, 'w+') do |f|  ## https://stackoverflow.com/questions/7911669/how-to-create-a-file-in-ruby
    f.puts("|*************************************************************************|")
    f.puts("|Report file of possible interactions between predicted co-expressed genes|")
    f.puts("|*************************************************************************|")
    f.puts("")
    #f.puts("Number of genes analysed (from the original txt): #{File.open(gene_list_txt, "r").each_line.count}")
    f.puts("Number of genes analysed (from the original txt): #{gene_txt.length()}")
    f.puts("#{all_networks.length()} networks involving genes in our original list have been detected. \n \n ")
    all_networks.each do |net|
      f.puts "Genes from the original list which interact with each other and form networks: #{net.genes_in_network.join(', ')} \n " # SHOWING THE GENES IN OUR ORIGINAL LIST THAT WILL INTERACT WITH EACH OTHER (FORMING NETWORKS)
      #f.puts net.network
      f.puts "Complete networks, including both the genes in the list (above) and other iterations with genes outside the list: \n "
      net.network.each do |key, value| # SHOWING FULL NETWORKS, WHICH WILL CONTAIN BOTH THE ORIGINAL GENES IN OUR LIST, AS WELL AS THOSE THAT ARE OUTSIDE THE LIST
          f.puts "#{key} has the following connections: #{value.join(', ')}"
      end
      f.puts '' # Printing an space between each Network found, to better visualize it in the terminal
      f.puts '_______________________________________________________________________________________________________'
      
    end
    
  end
end