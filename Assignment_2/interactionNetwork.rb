require 'rest-client'
require 'json'

class InteractionNetwork

    #code
    
    attr_accessor :network
    
    def initialize(params={})
      depth = params.fetch(:depth)
      @network = Hash.new
      gene_id = params.fetch(:gene_id)
      search_Interaction(gene_id, depth)
    end
    

    def search_Interaction(gene_id, depth)

      res = Utils.fetch ("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{gene_id}?species:3702?format=tab25")
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
          
          if g1.upcase != gene_id.upcase # Filter 2: Some gene positions in intact are swapped, we need to fix it
                g1, g2 = g2, g1
          end
          
          if g1 == g2 or g2 == "" #Filter 3: We must have into account that they could iterate with themselves, or with a non-arabidopsis gene
            next
          end
          
          if intact_miscore < 0.5#0.485 #Filter 4: When the miscore is lower than 0.485, this value was taken from the literature ##REVISAR
            next
          end
          # END OF FILTERS
          
          #Comprobar (Borrar)
          #puts "Pareja en-> #{gene_id} con #{g2}"
          
          
          #Comprobar si puedes guardarlo en el hash
          if @network.keys.include?(g2) #G2 ya existe en la network? 
            if @network[g2].include?(g1) #G2 ya tiene una relacion con g1?
                #This is a L00P!
                next
            end
          end
          
          
          #Guardarlo
          if @network.keys.include?(g1)
            #Si ya existe
            if @network[g1].include?(g2)
              next
            end
            @network[g1] << g2 #Append al array de g1
          else
            @network[g1] = [g2]
          end
          
          #Puedo seguir mirando? (DEPTH)
          if depth > 1
            #Seguir mirando (con un depth menos!)
            search_Interaction(g2, depth-1)
          end
          if depth == 0
            next
          end

        end
        # END OF FOR
        
      else
          puts "Error with: #{@gene_id}"
      end
        
    end
 
    
end