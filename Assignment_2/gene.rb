require 'rest-client'
require 'json'
require './utils.rb'

class Gene
    attr_accessor :gene_id
    attr_accessor :parejas

    
    def initialize(params={})
      @gene_id = params.fetch(:id)
      @parejas = Array.new
      #check_parejas()
    end



    def check_parejas()
      res = Utils.fetch ("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{@gene_id}?species:3702?format=tab25")
      
      if res
        lines = res.body.split("\n")
        #-------iteration of lines from ebi----------------

        lines.each do |i|
          line_splitted = i.split("\t")
          p1_locus = line_splitted[4]
          p2_locus = line_splitted[5]
          score = line_splitted[14]
          
          
          intact_miscore = score.sub(/intact-miscore:/, "").to_f #Extracting just the miscore number
          
          #Now we are going to apply some filters:
          p1_locus = p1_locus.match(/A[Tt]\d[Gg]\d{5}/).to_s.upcase #Filter 1:Applying a regular expression to match the Locus code
          if p1_locus == ""
            next
          else
            g1 = p1_locus
          end
         
          p2_locus = p2_locus.match(/A[Tt]\d[Gg]\d{5}/).to_s.upcase #Filter 1: Applying a regular expression to match the Locus code
          if p2_locus == ""
            next
          else
            g2 = p2_locus
          end
          #puts ""
          
          if p1_locus.upcase != @gene_id.upcase # Filter 2: Some gene positions in intact are swapped, we need to fix it
                p1_locus, p2_locus = p2_locus, p1_locus
          end
          
          if g1 == g2 or g2 == "" #Filter 3: We must have into account that they could iterate with themselves, or with a non-arabidopsis gene
            next
          end
          
          if intact_miscore < 0.485 #Filter 4: When the miscore is lower than 0.485, this value was taken from the literature ##REVISAR
            next
          end
          # END OF FILTERS
          
          #Comprobar (Borrar)
          puts "Pareja en-> #{@gene_id}"
          
          if @parejas.include? [[p1_locus, p2_locus, intact_miscore]]
            next
          end
          @parejas <<  [[p1_locus, p2_locus, intact_miscore]]
          
          #
    
        end
        # END OF FOR
        
      else
          puts "Error with: #{@gene_id}"
      end
        
    end



end

#URL EXAMPLE --> http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/AT2G47930?species:3702?format=tab25