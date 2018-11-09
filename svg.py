''' SVG.py - Construct/display SVG scenes.

The following code is a lightweight wrapper around SVG files. The metaphor
is to construct a scene, add objects to it, and then write it to a file
to display it. '''

class Scene:
    def __init__(self,name="svg",height=400,width=400):
        self.name = name
        self.items = []
        self.height = height
        self.width = width
    
    def add(self,item):
        self.items.append(item)
    
    def extend(self,items):
        self.items.extend(items)
    
    def strarray(self):
        var = ["<?xml version=\"1.0\"?>\n",
               "<svg height=\"%d\" width=\"%d\" >\n" % (self.height,self.width),
               " <g style=\"fill-opacity:1.0; stroke:black;\n",
               "  stroke-width:1;\">\n"] + \
               [item.strarray() for item in self.items] + \
               [" </g>\n</svg>\n"]
        return var
    
    def write_svg(self,filename=None):
        if filename:
            self.svgname = filename
        else:
            self.svgname = self.name + ".svg"
        with open(self.svgname,'w') as output:
            output.writelines(self.strarray())


class Line:
    def __init__(self,start,end):
        self.start = start  # xy tuple
        self.end = end  # xy tuple
    
    def strarray(self):
        return "  <line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" />\n" % \
                (self.start[0],self.start[1],self.end[0],self.end[1])


class Circle:
    def __init__(self,center,radius,color):
        self.center = center  # xy tuple
        self.radius = radius  # xy tuple
        self.color = color  # rgb tuple in range(0,256)
    
    def strarray(self):
        return "  <circle cx=\"%d\" cy=\"%d\" r=\"%d\"\n" % \
                (self.center[0],self.center[1],self.radius) + \
                "    style=\"fill:%s;\"  />\n" % colorstr(self.color)

class Rectangle:
    def __init__(self,origin,height,width,color):
        self.origin = origin
        self.height = height
        self.width = width
        self.color = color
    
    def strarray(self):
        return "  <rect x=\"%d\" y=\"%d\" height=\"%d\"\n" % \
                (self.origin[0],self.origin[1],self.height) + \
                "    width=\"%d\" style=\"fill:%s;\" />\n" % \
                (self.width,colorstr(self.color))

class Arc:
    def __init__(self,origin,end,radius,color,flip=False):
        self.origin = origin
        self.end = end
        self.radius = radius
        self.color = color
        self.flip = flip
    
    def strarray(self):
        return "  <path d=\"M %d,%d A %d,%d 0 0 %d %d,%d\"\n" % \
                (self.origin[0],self.origin[1],self.radius,self.radius,int(self.flip),self.end[0],self.end[1]) + \
                "    style=\"fill:none;stroke:%s;\"  />\n" % \
                colorstr(self.color)

class Text:
    def __init__(self,origin,text,size=24):
        self.origin = origin
        self.text = text
        self.size = size
    
    def strarray(self):
        return "  <text x=\"%d\" y=\"%d\" font-size=\"%d\">\n" % \
                (self.origin[0],self.origin[1],self.size) + \
                "   %s\n" % self.text + \
                "  </text>\n"

def colorstr(rgb): return "#%x%x%x" % (rgb[0]/16,rgb[1]/16,rgb[2]/16)

