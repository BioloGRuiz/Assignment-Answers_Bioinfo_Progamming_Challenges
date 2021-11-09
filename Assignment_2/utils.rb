require 'rest-client'

class Utils
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
  
  def self.read_txt(path)
  #How to read txt: https://www.rubyguides.com/2015/05/working-with-files-ruby/
  genes_list = File.read(path).split
  #How to put the whole array in caps: https://stackoverflow.com/questions/11402362/how-can-i-uppercase-each-element-of-an-array/11402608
  genes_list.map!(&:upcase)
  return genes_list
  end
end
