#!/usr/bin/env python3

#################################
# Import some necessary modules #
#################################

from optparse import OptionParser, OptionGroup

from Bio.Nexus import Trees, Nodes
from Bio.Graphics.GenomeDiagram._Colors import ColorTranslator
from Bio.GenBank import Scanner
from Bio.GenBank import _FeatureConsumer
from Bio.GenBank.utils import FeatureValueCleaner

from reportlab.lib.units import inch
from reportlab.lib import pagesizes
from reportlab.graphics.shapes import *
from reportlab.pdfgen.canvas import Canvas
from reportlab.graphics import renderPDF

################################
# Get the command line options #
################################


def main():
	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "Output Options")
	group.add_option("-o", "--output", action="store", dest="outputfile", help="output file name [default= %default]", type="string", metavar="FILE", default="test.pdf")
	group.add_option("-t", "--tree", action="store", dest="tree", help="tree file to align tab files to", default="")
	
	parser.add_option_group(group)
	
	return parser.parse_args()

  ##########################################################
  # Function to read an alignment whichever format it's in #
  ##########################################################


def tab_parser(handle, quiet=False):
	def Drawer_parse_tab_features(object, skip=False):
		features = []
		line = object.line
		while True:
			if not line:
				break
				raise ValueError("Premature end of line during features table")
			if line[:object.HEADER_WIDTH].rstrip() in object.SEQUENCE_HEADERS:
				if object.debug : print("Found start of sequence")
				break
			line = line.rstrip()
			if line == "//":
				raise ValueError("Premature end of features table, marker '//' found")
			if line in object.FEATURE_END_MARKERS:
				if object.debug : print("Found end of features")
				line = object.handle.readline()
				break
			if line[2:object.FEATURE_QUALIFIER_INDENT].strip() == "":
				print(line[2:object.FEATURE_QUALIFIER_INDENT].strip())
				raise ValueError("Expected a feature qualifier in line '%s'" % line)
			
			if skip:
				line = object.handle.readline()
				while line[:object.FEATURE_QUALIFIER_INDENT] == object.FEATURE_QUALIFIER_SPACER:
					line = object.handle.readline()
			else:
				#Build up a list of the lines making up this feature:
				feature_key = line[2:object.FEATURE_QUALIFIER_INDENT].strip()
				feature_lines = [line[object.FEATURE_QUALIFIER_INDENT:]]
				line = object.handle.readline()
				while line[:object.FEATURE_QUALIFIER_INDENT] == object.FEATURE_QUALIFIER_SPACER or line.rstrip() == "" : # cope with blank lines in the midst of a feature
					
					feature_lines.append(line[object.FEATURE_QUALIFIER_INDENT:].rstrip())
					line = object.handle.readline()
					if len(line)==0:
						break#EOF
				
				feature_lines.append('/seq="N"')
				sys.stdout.flush()
				features.append(object.parse_feature(feature_key, feature_lines))
		object.line = line
		
		return features

	
	def Drawer_feed(object, handle, consumer, do_features=True):
		if do_features:
			object._feed_feature_table(consumer, Drawer_parse_tab_features(object,skip=False))
		else:
			Drawer_parse_tab_features(object,skip=True) # ignore the data
		
		sequence_string="N"
		consumer.sequence(sequence_string)
		consumer.record_end("//")
		length=0
		for record in consumer.data.features:
			if record.location.nofuzzy_end>length:
				length=record.location.nofuzzy_end
		
		consumer.data.seq="N"*length
		
		return True
	
	myscanner=Scanner.InsdcScanner()
	myscanner.set_handle(handle)
	
	myscanner.line=myscanner.handle.readline()
	myscanner.FEATURE_QUALIFIER_INDENT=21
	myscanner.FEATURE_QUALIFIER_SPACER = "FT" + " " * (myscanner.FEATURE_QUALIFIER_INDENT-2)
	
	myscanner.debug=True
	
	consumer = _FeatureConsumer(use_fuzziness = 1, feature_cleaner = FeatureValueCleaner())
	
	Drawer_feed(myscanner, handle, consumer)
	
	return consumer.data

####################################################
# Function to round floats to n significant digits #
####################################################

def round_to_n(x, n):
    if n < 1:
        raise ValueError("number of significant digits must be >= 1")
    # Use %e format to get the n most significant digits, as a string.
    format = "%." + str(n-1) + "e"
    as_string = format % x
    if x>=10 or x<=-10:
    	return int(float(as_string))
    else:
	    return float(as_string)

##############################################################################################################
# Function to convert features with subfeatures (e.g. pseudogenes) to a list of locations of the subfeatures #
##############################################################################################################

def iterate_subfeatures(feature, locations):
	if len(feature.sub_features)>0:
		for subfeature in feature.sub_features:
			locations=iterate_subfeatures(subfeature, locations)
	else:
		locations.append((feature.location.start.position, feature.location.end.position))
	
	
	return locations


####################################################
# Function to get the pixel width of a text string #
####################################################

def get_text_width(font, size, text):
	c = Canvas(test,pagesize=pagesize)
	length= c.stringWidth(str(text),font,size)
	
	return length



#####################################################################################
# Function to add an embl file to multiple tracks split using the qualifiers option #
#####################################################################################

def add_ordered_embl_to_diagram(record, incfeatures=["CDS", "feature"], emblfile=True):
	incfeatures= [x.lower() for x in incfeatures]
	
	new_tracks={}
	
	print(len(record.features), "features found for", record.name)
	
	if len(record.seq)>500000:
		scale_largetick_interval=int(round((len(record.seq)/10),-5))
		scale_smalltick_interval=int(round((len(record.seq)/10),-5)/5)
	else:
		scale_largetick_interval=len(record.seq)
		scale_smalltick_interval=len(record.seq)/5
	
	for x, feature in enumerate(record.features):
		if feature.type.lower() not in incfeatures or feature.location.nofuzzy_end<0 or (feature.location.nofuzzy_start>-1 and -1!=-1):
			
			continue
		
		if "colour" in feature.qualifiers:
			colourline=feature.qualifiers["colour"][0]
		elif "color" in feature.qualifiers:
			colourline=feature.qualifiers["color"][0]
		else:
			colourline = "5"
		if len(colourline.split())==1:
			colour=translator.artemis_color(colourline)
		elif len(colourline.split())==3:
			colour=translator.int255_color((int(colourline.split()[0]),int(colourline.split()[1]),int(colourline.split()[2])))
		else:
			print("Can't understand colour code!")
			print(colourline)
			sys.exit()
		
		locations=[]
		

		locations.append((feature.location.start, feature.location.end))

		if "taxa" in feature.qualifiers:
			qualifiernames=feature.qualifiers["taxa"][0].replace(", "," ").split()
			
			for taxonname in qualifiernames:
				taxonname=taxonname.strip()
				if not taxonname in new_tracks:
					newtrack = Track()
					newtrack.name=taxonname
					new_tracks[taxonname]=newtrack
				
				arrows=0
				new_tracks[taxonname].add_feature(locations, fillcolour=colour, strokecolour=colour, arrows=arrows)
		
		else:
			if not record.name in new_tracks:
				newtrack = Track()
				newtrack.name=record.name
				new_tracks[record.name]=newtrack
			
			arrows=0
			new_tracks[record.name].add_feature(locations, fillcolour=colour, strokecolour=colour, arrows=arrows)
		
	
	if len(new_tracks)>1 and record.name in new_tracks:
		del new_tracks[record.name]
	return new_tracks

###################################################################################
# Function to add a tab file to multiple tracks split using the qualifiers option #
###################################################################################

def add_ordered_tab_to_diagram(filename):
	
	features={"":[]}
	
	featurename=""
	names_to_add_feature_to=[]
	
	try:
		record=tab_parser(open(filename,"r"))
	except IOError:
		print("Cannot find file", filename)
		sys.exit()
	record.name=filename
	new_tracks=add_ordered_embl_to_diagram(record, incfeatures=["i", "d", "li", "del", "snp", "misc_feature", "core", "cds", "insertion", "deletion", "recombination", "feature", "blastn_hit", "fasta_record", "variation"], emblfile=False)
	return new_tracks

def add_empty_track(existing_tracks, track_name):
	newtrack = Track()
	newtrack.name=track_name
	newtrack.beginning=0
	newtrack.track_height=1
	existing_tracks[track_name] = newtrack
	existing_tracks[track_name].add_feature(locations=[(0,0)], fillcolour=translator.artemis_color(2), strokecolour=translator.artemis_color(2), arrows=0)
	
	return existing_tracks
	
#############################
# Function to draw the tree #
#############################

def drawtree(treeObject, treeheight, treewidth, xoffset, yoffset, name_offset=5):
	
	
	def get_max_branch_depth():
		
		terminals=treeObject.get_terminals()
		maxbrlen=0.0
		for terminal in terminals:
			if treeObject.sum_branchlength(node=terminal)>maxbrlen:
				maxbrlen=treeObject.sum_branchlength(node=terminal)
		
		return maxbrlen
	
	def draw_scale():
		
		if vertical_scaling_factor<5:
			linewidth=0.5
		else:
			linewidth=1.0
		branchlength=round_to_n(max_branch_depth/10, 2)*horizontal_scaling_factor
		horizontalpos=xoffset+round_to_n(max_branch_depth/10, 2)*horizontal_scaling_factor
		vertpos=treebase-fontsize
		scalestring = str(round_to_n(max_branch_depth/10, 2))
		scalefontsize=fontsize
		if scalefontsize<6:
			scalefontsize=6
		d.add(Line(horizontalpos, vertpos, horizontalpos+branchlength, vertpos, strokeWidth=linewidth))
		d.add(String(horizontalpos+(float(branchlength)/2), vertpos-(scalefontsize+1), scalestring, textAnchor='middle', fontSize=scalefontsize, fontName='Helvetica'))

	
	
	def get_node_vertical_positions():
		
		def get_node_vertical_position(node):
			
			for daughter in treeObject.node(node).succ:
				get_node_vertical_position(daughter)

			
			if not treeObject.is_terminal(node):
				daughters=treeObject.node(node).succ
				if treeObject.node(node).data.comment==None:
					treeObject.node(node).data.comment={}
				treeObject.node(node).data.comment["vertpos"]=float(treeObject.node(daughters[0]).data.comment["vertpos"]+treeObject.node(daughters[-1]).data.comment["vertpos"])/2
		
		
		node=treeObject.root
		get_node_vertical_position(node)
	
	
	def drawbranch(node,horizontalpos):
		
		vertpos=treeObject.node(node).data.comment["vertpos"]+yoffset
		
		horizontalpos+=xoffset
		
		branchlength=treeObject.node(node).data.branchlength*horizontal_scaling_factor
		
		if vertical_scaling_factor<5:
			linewidth=0.5
		else:
			linewidth=1.0
		
		if treeObject.node(node).data.comment and "branch_colour" in treeObject.node(node).data.comment:
			r,g,b=treeObject.node(node).data.comment["branch_colour"]
			branch_colour=colors.Color(float(r)/255,float(g)/255,float(b)/255)
		else:
			branch_colour=colors.black
		if branchlength<linewidth:
			branchlength=linewidth
		d.add(Line(horizontalpos-(linewidth/2), vertpos, (horizontalpos-(linewidth/2))+branchlength, vertpos, strokeWidth=linewidth, strokeColor=branch_colour))
		
		
		if node!=treeObject.root:
			
			parentnode=treeObject.node(node).prev
			sisters=treeObject.node(parentnode).succ
			parentvertpos=treeObject.node(parentnode).data.comment["vertpos"]+yoffset
			d.add(Line(horizontalpos, vertpos, horizontalpos, parentvertpos, strokeWidth=linewidth, strokeColor=branch_colour))
		
		
		if treeObject.is_terminal(node):
			
			if treeObject.node(node).data.comment and "name_colour" in treeObject.node(node).data.comment:
				name_colours=[]
				for x in range(0,len(treeObject.node(node).data.comment["name_colour"])):
					r,g,b= treeObject.node(node).data.comment["name_colour"][x]
					name_colours.append(colors.Color(float(r)/255,float(g)/255,float(b)/255))
			else:
				name_colours=[colors.black]
			
			gubbins_length=0.0
			
			colpos=0
			namewidth=get_text_width('Helvetica', fontsize, treeObject.node(node).data.taxon)+name_offset
			gubbins_length += namewidth
			colpos=1
			
			for x in range(colpos,len(name_colours)):
				gubbins_length += block_length
				if x!=0:
					gubbins_length += vertical_scaling_factor
			
			#Add the taxon names
			d.add(String(treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2), vertpos-(fontsize/3), treeObject.node(node).data.taxon, textAnchor='start', fontSize=fontsize, fillColor=name_colours[0], fontName='Helvetica'))
			block_xpos=treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2)+namewidth
			
			# draw dashed lines
			d.add(Line(horizontalpos+branchlength, vertpos, treewidth+xoffset+(max_name_width-gubbins_length), vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colours[0]))

	
	def recurse_subtree(node, horizontalpos):
		
		daughters=treeObject.node(node).succ
		
		daughterhorizontalpos=horizontalpos+(treeObject.node(node).data.branchlength*horizontal_scaling_factor)
		drawbranch(node,horizontalpos)
		for daughter in daughters:
			recurse_subtree(daughter,daughterhorizontalpos)
		

	def get_max_name_width(name_offset, fontsize):
		max_width=0.0
		for taxon in treeObject.get_terminals():
			curwidth= get_text_width("Helvetica", fontsize, treeObject.node(taxon).data.taxon)
			if curwidth>max_width:
				max_width=curwidth
		
		return max_width

	fontsize=vertical_scaling_factor
	if fontsize>12:
		fontsize=12
	
	while get_max_name_width(name_offset, fontsize)+name_offset>treewidth/3:
		fontsize-=0.2
	max_name_width=get_max_name_width(name_offset, fontsize)+name_offset
	colblockstart=1
	
	block_length=0
	
	treewidth-=(max_name_width+(fontsize/2))
	max_branch_depth=get_max_branch_depth()
	horizontal_scaling_factor=float(treewidth)/max_branch_depth
	get_node_vertical_positions()
	recurse_subtree(treeObject.root, 0)
	treebase=treeObject.node(treeObject.get_terminals()[-1]).data.comment["vertpos"]+yoffset
	
	draw_scale()
	
	return



#################
# Drawing class #
#################

class Figure:
	def __init__(self, beginning, end):
		self.begnining=0
		self.end=-1


###############
# Track class #
###############


class Track:
	def __init__(self, track_position=[-1,-1], track_height=0, track_length=0, track_draw_proportion=0.75, scale=False, tick_marks=True, tick_mark_number=5, tick_mark_labels=True, minor_tick_marks=True, minor_tick_mark_number=3, features=[], beginning=0, end=-1):
		
		self.track_position=track_position#horizontal and vertical position of centre of track
		self.track_height=track_height#height of space allocated for track
		self.track_length=track_length
		self.track_draw_proportion=track_draw_proportion#proportion of the track that should be used for drawing
		self.scale=scale
		self.scale_position="middle"
		self.tick_marks=tick_marks
		self.tick_mark_number=tick_mark_number
		self.tick_mark_labels=tick_mark_labels
		self.tick_mark_label_font="Helvetica"
		self.tick_mark_label_size=8
		self.tick_mark_label_angle=45
		self.minor_tick_marks=minor_tick_marks
		self.minor_tick_mark_number=minor_tick_mark_number
		self.features=features[:]
		self.scaled_features=features[:]
		self.draw_feature_labels=False
		self.feature_label_size=8
		self.feature_label_angle=0
		self.feature_label_font="Helvetica"
		self.greytrack=False
		self.grey_track_colour=colors.Color(0.25,0.25,0.25)
		self.grey_track_opacity_percent=10
		self.max_feature_length=-1
		self.beginning=0
		self.end=-1
		self.track_number=-1
		self.plots=[]
		self.fragments=1
		self.name=""
		self.show_name=False
		self.name_font="Helvetica"
		self.name_size=10
		self.name_length=0
		self.is_key=False
		self.key_data=[]
	
	
	def get_max_feature_length(self):
		max_feature_length=0
		for feature in self.features:
			for location in feature.feature_locations:
				if location[0]>max_feature_length:
					max_feature_length=location[0]
				if location[1]>max_feature_length:
					max_feature_length=location[1]
		return max_feature_length
	
	
	def scale_feature_positions(self):
		
		self.scaled_features=[]
		
		if self.end!=-1:
			length=float(self.end-self.beginning)
		else:
			length=float(self.max_feature_length-self.beginning)
		
		for feature in self.features:
			
			newfeature=Feature()
			newfeature.fillcolour=feature.fillcolour
			newfeature.strokecolour=feature.strokecolour
			newfeature.strokeweight=feature.strokeweight
			newfeature.strand=feature.strand
			newfeature.label=feature.label
			newfeature.arrows=feature.arrows
			scaledlocations=[]
			for location in feature.feature_locations:
				start=location[0]
				finish=location[1]
				if self.beginning!=0:
					if start<self.beginning and finish>self.beginning:
						start=self.beginning
				if self.end!=-1:
					if start<self.end and finish>self.end:
						finish=self.end
				start-=self.beginning
				finish-=self.beginning
				
				scaledlocations.append(((float(start)/length)*self.track_length,(float(finish)/length)*self.track_length))
			
			newfeature.feature_locations=scaledlocations
			self.scaled_features.append(newfeature)
		
	
	def draw_features(self):
		
		if self.max_feature_length==-1:
			return
		
		else:
			self.scale_feature_positions()
		
		featuresort=[]
		for x, feature in enumerate(self.scaled_features):
			featuresort.append([feature.feature_locations[0][0], x])
		joins=[]
		
		
		
		for featurenum in featuresort[::-1]:
			feature=self.scaled_features[featurenum[1]]
			#if the feature is white, outline it in black so we can see it
			if feature.strokecolour==colors.Color(1,1,1,1):
				feature.strokecolour=colors.Color(0,0,0,1)
			
			subfeaturesort=[]
			for x, subfeature in enumerate(feature.feature_locations):
				subfeaturesort.append([subfeature[0], x])
			subfeaturesort.sort()
			subfeature_locations=[]
			for subfeaturenum in subfeaturesort:
				subfeature_locations.append(feature.feature_locations[subfeaturenum[1]])

			for x, location in enumerate(subfeature_locations):
				
				
				if (location[0]>0 and location[0]<=self.track_length) or (location[1]>0 and location[1]<=self.track_length):
					
					y=self.track_position[1]-((float(self.track_height)/4)*self.track_draw_proportion)
					height=(float(self.track_height)*self.track_draw_proportion)/2
					y1=self.track_position[1]
					y2=self.track_position[1]+((float(self.track_height)/8)*self.track_draw_proportion)
					
					
					if feature.arrows==0:
						d.add(Rect(self.track_position[0]+location[0], y, location[1]-location[0], height, fillColor=feature.fillcolour, strokeColor=feature.strokecolour, strokeWidth=feature.strokeweight))
					
					if len(subfeature_locations)>x+1 and subfeature_locations[x+1][0]<=self.track_length:
						if subfeature_locations[x+1][0]<location[1]:
							joinheight=y1
						elif y2>y1:
							if (y2-y1)>(float(subfeature_locations[x+1][0]-location[1])/2):
								joinheight=y1+(float(subfeature_locations[x+1][0]-location[1])/2)
							else:
								joinheight=y2
						else:
							if (y1-y2)>(float(subfeature_locations[x+1][0]-location[1])/2):
								joinheight=y1-(float(subfeature_locations[x+1][0]-location[1])/2)
							else:
								joinheight=y2
						
						joins.append(Line(self.track_position[0]+location[1], y1, self.track_position[0]+location[1]+(float(subfeature_locations[x+1][0]-location[1])/2), joinheight, strokeDashArray=[0.5, 1], strokeWidth=0.5))
						joins.append(Line(self.track_position[0]+((location[1]+subfeature_locations[x+1][0])/2), joinheight, self.track_position[0]+location[1]+(float(subfeature_locations[x+1][0]-location[1])), y1, strokeDashArray=[0.5, 1], strokeWidth=0.5))
				
		
		for join in joins:
			d.add(join)
		
		self.scaled_features=[]

	
	def draw_track(self):
		self.draw_features()

	
	def add_feature(self,locations=[(-1,-1)], fillcolour=colors.white, strokecolour=colors.black, strokeweight=0, label="", strand=0, arrows=0):
		
		newfeature=Feature()
		
		feature_locations=[]
		for location in locations:
			if location[0]>location[1]:
				feature_locations.append((location[1],location[0]))
			else:
				feature_locations.append((location[0],location[1]))
		
		newfeature.feature_locations=feature_locations
		newfeature.fillcolour=fillcolour
		newfeature.strokecolour=strokecolour
		newfeature.strokeweight=strokeweight
		newfeature.strand=strand
		newfeature.label=label
		newfeature.arrows=arrows
		
		
		self.features.append(newfeature)
		
	
	def sort_features_by_length(self):
		featurelist=[]
		ordered_features=[]
		for x, feature in enumerate(self.features):
			featurelist.append([feature.feature_locations[-1][1]-feature.feature_locations[0][0], x])
		featurelist.sort()
		#featurelist.reverse()
		
		for feature in featurelist:
			ordered_features.append(self.features[feature[1]])
		self.features=ordered_features[:]


#################
# Feature class #
#################


class Feature:
	def __init__(self):
		
		self.feature_locations=[(-1,-1)]
		self.strand=0
		self.arrows=0
		self.label=""
		self.fillcolour=colors.blue
		self.strokecolour=colors.black
		self.strokeweight=0


################
# Main program #
################

if __name__ == "__main__":

	(options, args) = main()
	pagesize=pagesizes.A4
	height, width = pagesize
  
	if len(args)==0:
		print("Found nothing to draw")
		sys.exit()
	
	d = Drawing(width, height)
	margin=0.5*inch
	metadatanames={}
	namecolours={}
	colour_dict=[]
	my_tracks={}
	
	#create translator object for translating artemis colours to GenomeDiagram colours
	translator = ColorTranslator()
	track_count=0
	tree_positions=[]
	track_names={}
	input_order=[]
	
	
	for arg in args[::-1]:
		if arg.lower() in ["tree", "list"]:
			input_order.append(arg.lower())
			continue
		if arg.split('.')[-1].lower() in ["plot", "hist", "heat", "bar", "line", "graph", "area","embl", "gb", "tab", "bam", "fas", "fasta", "mfa", "dna", "fst", "phylip", "phy", "nexus", "nxs"]:
			
			if arg.split('.')[-1].lower() =="embl":
				 new_tracks=add_ordered_tab_to_diagram(arg)
				 
				 for track in new_tracks:
				 	newtrack=new_tracks[track]
				 	newtrack.beginning=0
				 	newtrack.name=new_tracks[track].name
				 	name=newtrack.name
				 	x=1
				 	while name in my_tracks:
				 		name=newtrack.name+"_"+str(x)
				 		x+=1
				 	if not newtrack.name in track_names:
				 		track_names[newtrack.name]=[]
				 	input_order.append(name)
				 	track_names[newtrack.name].append(name)
				 	
				 	track_count+=1
				 	newtrack.track_height=1
				 	my_tracks[name]=newtrack
	
	treenames=[]
	tree_name_to_node={}
	listnames=[]
	
	if options.tree!="":
		if not os.path.isfile(options.tree):
			print("Cannot find file:", options.tree)
			options.tree=""
		else:
			treestring=open(options.tree,"rU").read().strip()
			tree=Trees.Tree(treestring, rooted=True)
			
			tree.root
			
			treeterminals=tree.get_terminals()
			totalbr=0.0
			
			
			for terminal_node in treeterminals:
				terminal=tree.node(terminal_node).data.taxon
				treenames.append(terminal)
				if not terminal in track_names:
					track_count+=1
				tree_name_to_node[terminal]=terminal_node
				tree.node(terminal_node).data.comment={}
				tree.node(terminal_node).data.comment["name_colour"]=[(0,0,0)]

	
	
	#from this we can work out a constant for the height of a track which takes into account the height of the page and margin sizes
	
	vertical_scaling_factor=float(height-(margin*2))/(track_count)
	
	#to make sure names can be printed in the space of a track, we can scale the name to the same size as the vertical scaling factor, but limit it to 12pt so it doesn't get crazily big
	
	name_font_size=vertical_scaling_factor
	if name_font_size>12:
		name_font_size=12
	
	left_proportion=0.3

	treetrack=0
	output_order=treenames[::-1]
	
	
	for name in input_order[::-1]:
		if not name in treenames:
			output_order.append(name)
	
	track_number=0
	
	for track in output_order:
		if(track not in my_tracks):
			my_tracks = add_empty_track(my_tracks, track)
		
		track_height=my_tracks[track].track_height
		
		my_tracks[track].track_draw_proportion=0.8
		my_tracks[track].track_height=track_height*vertical_scaling_factor
		if left_proportion==1:
			my_tracks[track].track_length=(width-margin)-((width-(margin*2))*0.2+margin)
			my_tracks[track].track_position=[(width-(margin*2))*0.2+margin, margin+((track_number)*vertical_scaling_factor)+float((my_tracks[track].track_height)/2)]
		
		else:
			my_tracks[track].track_length=(width-margin)-((width-(margin*2))*left_proportion+margin)
			my_tracks[track].track_position=[(width-(margin*2))*left_proportion+margin, margin+((track_number)*vertical_scaling_factor)+float((my_tracks[track].track_height)/2)]
		
		my_tracks[track].track_number=track_number
		if track in treenames:
			tree.node(tree_name_to_node[track]).data.comment["vertpos"]=margin+((track_number)*vertical_scaling_factor)+float((my_tracks[track].track_height)/2)
			my_tracks[track].grey_track_colour=colors.Color(0,0,0)
		
		track_number+=track_height


	
	#find the maximum feature endpoint to scale by
	max_feature_length=0
	
	for track in my_tracks:
		max_track_feature_length=my_tracks[track].get_max_feature_length()
		if max_track_feature_length>max_feature_length:
			max_feature_length=max_track_feature_length
		for plot in my_tracks[track].plots:
			for data in plot.xdata:
				if data[-1]>max_feature_length:
					max_feature_length=data[-1]
	
	#tell each track what the max feature length is
	for track in my_tracks:
		if my_tracks[track].max_feature_length<max_feature_length:
				my_tracks[track].max_feature_length=max_feature_length


	beginning=0
	end=max_feature_length
	
	for track in output_order:
		
		if not track in my_tracks or (my_tracks[track].is_key and fragment!=1) or my_tracks[track].track_length==0:
			continue
		
		my_tracks[track].beginning=beginning
		my_tracks[track].end=end
		
		my_tracks[track].track_position[1]=margin+(((my_tracks[track].track_number)*vertical_scaling_factor)+(my_tracks[track].track_height)/2)
		
		my_tracks[track].sort_features_by_length()
		my_tracks[track].draw_track()
	
	if options.tree!="":
		drawtree(tree, height-(margin*2), (width-(margin*2))*left_proportion, margin, 0, 5)
	
	renderPDF.drawToFile(d, options.outputfile)

class DrawerError(Exception):
  def __init__(self, value):
    		self.value = value
  def __str__(self):
    		return repr(self.value)

