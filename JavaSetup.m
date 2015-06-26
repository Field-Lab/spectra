% Work in progress utilitarian file
% Contents can change and should be user defined

%% Java import setup
javaaddpath('C:/Users/Vincent/Documents/EJGroup/Java/vision7/bin/');
javaaddpath('C:/Users/Vincent/Documents/EJGroup/Java/vision7/Vision.jar','-end');
import edu.ucsc.neurobiology.vision.*
import edu.ucsc.neurobiology.vision.anf.*
import edu.ucsc.neurobiology.vision.math.*
import edu.ucsc.neurobiology.vision.calculations.*
import edu.ucsc.neurobiology.vision.io.*
import edu.ucsc.neurobiology.vision.electrodeMap.*
import java.io.*
import java.util.*

% Set up test path
dataPath = 'X:\EJGroup_data\Data\2005-04-26-0\data002';
saveFolder = 'X:\EJGroup_data\TestOut\2005-04-26-0\data002';

% Set up vision config XML file
config = edu.ucsc.neurobiology.vision.Config('C:\Users\Vincent\Documents\EJGroup\Java\vision7\config.xml');

%% Set up a HashMap of parameters for SpikeFinding
sigmaFileName = [saveFolder,'\data002.noise'];
spikeFileName = [saveFolder,'\data002.spikes'];
parameters = spikeFindingSetup([dataPath,'(0-10)'],saveFolder,sigmaFileName,config);