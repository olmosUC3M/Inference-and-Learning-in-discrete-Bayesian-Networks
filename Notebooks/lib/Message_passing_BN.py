## Message passing over a discrete BN ##
## Library created by Pablo Martínez Olmos, University Carlos III Madrid ##
## olmos@tsc.uc3m.es ##
## Last modification 15/11/2016 ##

import numpy as np

## Messages are stored in the logaritmic domain ##
## Global constants (to control numerical issues)

inf_log=100		#To impose hard constraints (i.e. an observed variable)
constant_log=50		#Used to improve stability in the Check Node (CN) operation


## Function definitions


def create_var_node(ID,cardinality,neighbor_order,observed_value_index=-1):

    # Variable Nodes are defined by a dictionary with several fields	
    var_node={}
    var_node['ID']=ID
    var_node['node_type']=0	#type 0 refers to variable node, 1o to check nodes. 
    var_node['cardinality']=cardinality		#Num. of possible values the RV can take
    var_node['neighbor_order']=np.array(neighbor_order)	#Ordered array of the neighbor's IDs (neighbors are CNs!)
    var_node['input_msgs']=[]	#List to store input messages 
    var_node['observed']=observed_value_index	#-1 if the variable is not observed
    var_node['inner_factor']=np.zeros([cardinality,1])	#Internal vector used to imposed hard messages when variable is observed
    
    #If variable is observed, then the inner_factor vector is log[0 0 ... 0 1 0 ...]
    if(observed_value_index!=-1):
        var_node['inner_factor']-=inf_log
        var_node['inner_factor'][observed_value_index]=inf_log
          
    #Initialize input msgs by filling with zeros
    
    for index,f in enumerate(var_node['neighbor_order']):
        var_node['input_msgs'].append(0)
    
    return var_node

def create_message(input_node,output_node,table):

    #Messages are defined by a dictionary with three keys: input node (sender node), output_node (receiver node), and table of values
   
    message={}
    message['input_node']=input_node
    message['output_node']=output_node
    message['table']=table
    
    return message
  

def create_factor_node(ID,neighbors,CPD):
    
    
    # Check Nodes are defined by a dictionary with several fields	
    
    factor_node={}
    factor_node['ID']=ID
    factor_node['node_type']=1
    factor_node['input_msgs']=[]

    CPD=np.array(CPD)
    CPD=CPD.reshape(CPD.shape[0],)	#Just to make sure that CPD is a np. array vector of dim. (n,)
    factor_node['CPD']=np.array(CPD)	#CPD table associated to the factor

    factor_node['CPD_order']=np.zeros([len(neighbors),1]).astype(int) #Ordered array of the neighbor's IDs (neighbors are CNs!)
    factor_node['cardinalities']=np.zeros([len(neighbors),1]).astype(int) #Cardinalities of the neighbors
    

    #Initialize input msgs, CPD_order & cardinalities
    #Note that creating factor nodes requires variable nodes to be created first, as CN input messages 
    #are initialized already to the inner_factor field of every neighbor variable node

    for index,node in enumerate(neighbors):
        card=node['cardinality']
        factor_node['input_msgs'].append(
            create_message(input_node=node,output_node=factor_node,table=node['inner_factor']))
        factor_node['cardinalities'][index]=card
        factor_node['CPD_order'][index]=node['ID']
        
    return factor_node


def initialize_variable(var_node,observed_value_index=-1):

    #After running message passing, variable nodes store the incoming messages for future calculations
    #If we want to run again message passing in the same graph, we have to re-initialize both
    #variable nodes and check nodes.

    var_node['inner_factor']=np.zeros([var_node['cardinality'],1])
    var_node['observed']=observed_value_index

    if(observed_value_index!=-1):
        var_node['inner_factor']-=inf_log
        var_node['inner_factor'][observed_value_index]=inf_log
    
def initialize_factor_msgs(factor_node,neighbors):

    #After running message passing, variable nodes store the incoming messages for future calculations
    #If we want to run again message passing in the same graph, we have to re-initialize both
    #variable nodes and check nodes.

    factor_node['input_msgs']=[]    
    
    for index,node in enumerate(neighbors):
        factor_node['input_msgs'].append(
            create_message(input_node=node,output_node=factor_node,table=node['inner_factor']))
        

    #The next two routines are used to encode and decode positions to store CPD values in a
    #vector form. We use a tree-encoding determined by the order of variables and their cardinalities
    #See First Example Message Passing.ipynb for an illustration 

def CPD_position_to_variable_index(position,v_card,CPD_size):

    #We use this function to find the encoding for each position of a CPD table
    #of CPD_size positions, where the cardinalities of the variables (in order) are given in v_card 
    #This function returns the index value of each variable  

    v_card=np.array(v_card)   #To make sure we have a np.array

    var_index=np.zeros([v_card.shape[0],1]).astype(int)
    
    remaining=CPD_size
    for i,card in enumerate(v_card):
        remaining=remaining//card
        index_i=position//remaining
        position=position-index_i*(remaining)
        var_index[i]=index_i
            
    return var_index

def variable_index_to_CPD_position(var_index,v_card,CPD_size):

    #This function returns the encoded CPD position for a given configuration of the variables. 
    #The CPD table is of size CPD_size, the cardinalities of the variables (in order) are given in v_card
    #and the value indexes (in order) of the variables are given in var_index

    var_index=np.array(var_index)
    v_card=np.array(v_card)
    
    position=0
    offset=CPD_size
    for i,card in enumerate(v_card):
        offset=offset//card
        position+=var_index[i]*offset
    return position


def update_var_to_factor(var_node):
    
    #Routine to update the output messages of a variable node (var_node)

    prod_table=np.zeros([var_node['cardinality'],1])

    #We first multiply all the input messages (sums in the log domain)
    for msg in var_node['input_msgs']:
        prod_table+=msg['table']

    #We also take into account the inner_factor of the variable_node. In
    #case it is observed, the output messages have to be consistent with the observation
    prod_table+=var_node['inner_factor']
    

    #For every output message, we have to substract from prod_table the message received  
    #through the corresponding edge

    for msg in var_node['input_msgs']:

        if(var_node['observed']==-1):
            reply_table=prod_table-msg['table']
        else:
            reply_table=np.ones([var_node['cardinality'],1])*(-inf_log)
            reply_table[var_node['observed']]=inf_log

        #We limit the absolute value of the messages, to exp(inf_log)

        reply_table[reply_table>inf_log]=inf_log
        reply_table[reply_table<-inf_log]=-inf_log

	#The ouput message is stored in the corresponding neighbor
        factor_rx=msg['input_node']
        reply_msg=create_message(input_node=var_node,output_node=factor_rx,table=reply_table)
        
        #Short foor loop to save messages in factor_node in the corresponding order
        for index,v in enumerate(factor_rx['CPD_order']):
             if(v==var_node['ID']):
                factor_rx['input_msgs'][index]=reply_msg
                break



def compute_var_marginal(var_node):
    
    #Routine to compute the marginal pmf of a variable node (var_node)
    #Simply the product of all incoming msgs times the inner_factor

    marg_table=np.zeros([var_node['cardinality'],1])

    for msg in var_node['input_msgs']:
        marg_table+=msg['table']

    marg_table+=var_node['inner_factor']
                    
    marg_table=np.exp(marg_table)
    marg_table/=sum(marg_table)

        
    return marg_table


def update_factor_to_var(factor_node):

   #Routine to update the output messages of a check node (var_node)
   #This is the most complicated in the library, as it involves marginalization
   #over each argument of the CPD function times the product of incoming messgaes

 
    output_tables=[]
    
    #Output message tables initialization    
    for card in factor_node['cardinalities']:
        output_tables.append(np.zeros([card,1]))
    

    #With a single loop we go only once through every element of the CPD table
    #It is multiplied accordingly to input messages and the resulting terms are
    #added to the corresponding output tables
    
    for CPD_entry,CPD_val in enumerate(factor_node['CPD']):
        
        values=CPD_position_to_variable_index(
            position=CPD_entry,v_card=factor_node['cardinalities'],CPD_size=factor_node['CPD'].shape[0])
        
	#The CPD value is multiplied by all incoming input messages but one, 
	#and the result is added to the ouput table

	#Since we have to marginalize, not all operations can be done in the log domain
	#To avoid numerical inestabilities when performing the operations, we substract a large exponent (constant log)
	#which is sum at the very end, when we move back to the log domain

        for index in range(factor_node['cardinalities'].shape[0]):
            
            aux=CPD_val
            for index2 in range(factor_node['cardinalities'].shape[0]):
                if(index2!=index):
                    aux*=np.exp(factor_node['input_msgs'][index2]['table'][values[index2]]-constant_log)
            output_tables[index][values[index]]+=aux
    
    #Once the output tables have been computed, we create the output messages and store them in 
    #the corresponding variable nodes
    
    for index,msg in enumerate(factor_node['input_msgs']):
        
        output=output_tables[index]
        output=np.log(output)+constant_log
        output[output>inf_log]=inf_log
        output[output<-inf_log]=-inf_log
        
        var_rx=msg['input_node']
        reply_msg=create_message(input_node=factor_node,output_node=var_rx,table=output)
        
        #Short foor loop to save messages in factor_node in the corresponding order
        for index2,f in enumerate(var_rx['neighbor_order']):
            if(f==factor_node['ID']):
                var_rx['input_msgs'][index2]=reply_msg
                break
              
       
    
def create_joint_node(ID,node_members,neighbor_order,observed_values_indexes=-1):

    #Routine to define a joint variable node. This is useful to eliminate cycles in
    #the factor graph and perform exact inference.

    #Note a routine to create a joint factor node that uses joint variable nodes
    #is not provided. The corresponding CPD of such factor nodes has to be computed
    #first and then create the joint node with the function create_factor_node

    #We do not consider the case that the joint variable node is partially observed 
    #(e.g. one of the joined variable nodes is observed). We only consider the case
    #where the joint node is completely observed.

    #See Second Example Message Passing.ipynb for an example of how to define and 
    #manage joint variable nodes.

    var_node={}
    var_node['ID']=ID
    var_node['node_type']=0
    var_node['input_msgs']=[]
    var_node['observed']=-1
    var_node['neighbor_order']=np.array(neighbor_order)
    
    card=1

    #Cardinality of joint node is the product of cardinalities
    for member in node_members:
        card*=member['cardinality']
        
    var_node['cardinality']=card
    
    var_node['inner_factor']=np.zeros([card,1])
    
    if(observed_values_indexes!=-1):
        var_node['observed']=variable_index_to_CPD_position(observed_values_indexes,var_node['values'],card)
        var_node['inner_factor']-=inf_log
        var_node['inner_factor'][var_node['observed']]=inf_log
        
    #Initialize input msgs
    
    for index,f in enumerate(var_node['neighbor_order']):
        var_node['input_msgs'].append(0)    
    
    return var_node  
    
    
