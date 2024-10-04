#include "CParticleCycling.h"

//use the necessary input text file  : edit:abhilash
const string CParticleCycling::PREFS_FILE = "ParticleCycleInput.txt";
const string CParticleCycling::MAX_SCALE = "Initial particle Scale";

// create new property to store the cycle number : edit : Abhi
// /!\ these properties must also be declared in the header file CParticleCycling.h
// indexed and initialised in "CParticleCycling::getDetailsForProperty" and "CParticleCycling::starting" 
// 
// Working:
// a value of m_maxcycles is taken from the input prefs file, indicating the number of cycles to be run.
// the property index is extracted and stored in variable iCycle_Number the value of it is initialised to zero
// for the calculations the value is extracted by indexing to the property number and stored in variable CycleIndex
// CycleIndex is used to change the cycling conditions.(see calculations)
// once the end of first cycle is reached and the cycleindexdelta variable is incremented

const string CParticleCycling::Cycle_Number = "Cycle Number";


CParticleCycling::CParticleCycling()
{
    ;
}

CParticleCycling::~CParticleCycling()
{
    ;
}

bool CParticleCycling::isThreadSafe()
{
    // thread safe
    return true;
}

bool CParticleCycling::usesCustomProperties()
{
    // Uses custom properties
    return true;
}

void CParticleCycling::getPreferenceFileName(char prefFileName[NApi::FILE_PATH_MAX_LENGTH])
{
	strcpy(prefFileName, PREFS_FILE.c_str());
}

bool CParticleCycling::setup(NApiCore::IApiManager_1_0 & apiManager, const char prefFile[], char customMsg[NApi::ERROR_MSG_MAX_LENGTH])
{
	ifstream myinputfile(prefFile);

	if (!myinputfile)
	{
		strcpy(customMsg, "There is no input file");
		return false;
	}
	string thisLine;
	int dataNumber = 0;
	m_particleNames.clear();
	int nbOfParticleType;

	while (!myinputfile.eof())
	{
		// Skip lines with white space lines
		bool dataLine = false;
		//clears the previously stored line if any
		thisLine.clear(); 
		getline(myinputfile, thisLine);
		int pos = thisLine.find_first_not_of("\t\v\r\n''");

		// if a line has characters on it, establish if is is a comment
		if (pos != -1)
		{
			pos = thisLine.find("#");
			if (pos == -1) { dataLine = true; }
		}


		if (dataLine)
		{
			switch (dataNumber)
			{
			default:
				break;


			case(0) :
			{
				istringstream ss(thisLine);
				ss >> nbOfParticleType;
				dataNumber++;
				break;
			}

			case(1):
			{
				// store the property max number of cycles : edit : Abhi
				istringstream ss(thisLine);
				ss >> m_maxcycles;
				dataNumber++;
				break;
			}

			case(2) :
				m_particleNames.push_back(thisLine);
				for (int i = 0; i < nbOfParticleType - 1; ++i)
				{
					string temp;
					myinputfile >> temp;
					m_particleNames.push_back(temp);
					
				}
				dataNumber++;
				break;
			case(3) :
			{
				double temp;
				istringstream ss(thisLine);
				ss >> temp;
				m_expansionType.push_back(temp);
				for (int i = 0; i < nbOfParticleType - 1; ++i)
				{
					myinputfile >> temp;
					m_expansionType.push_back(temp);
				}
				dataNumber++;
				break;
			}
			case(4) :
			{
				double temp;
				istringstream ss(thisLine);
				ss >> temp;
				m_expansionFinalScale.push_back(temp);
				for (int i = 0; i < nbOfParticleType - 1; ++i)
				{
					myinputfile >> temp;
					m_expansionFinalScale.push_back(temp);
				}
				dataNumber++;
				break;
			}
			case(5) :
			{
				double temp;
				istringstream ss(thisLine);
				ss >> temp;
				m_expansionStartTime.push_back(temp);
				for (int i = 0; i < nbOfParticleType - 1; ++i)
				{
					myinputfile >> temp;
					m_expansionStartTime.push_back(temp);
				}
				dataNumber++;
				break;
			}
			//DwellTime accept the time to swell and deswell : edit:abhilash
			case(6):
			{
				double temp;
				istringstream ss(thisLine);
				ss >> temp;
				m_DwellTime.push_back(temp);
				for (int i = 0; i < nbOfParticleType - 1; ++i)
				{
					myinputfile >> temp;
					m_DwellTime.push_back(temp);
				}
				dataNumber++;
				break;
			}

			}
		}
	}
	
	return true;
}
bool CParticleCycling::starting(NApiCore::IApiManager_1_0& apiManager, int numThreads)
{
    m_particleMngr = static_cast<NApiCore::IParticleManagerApi_1_0*>(apiManager.getApi(eParticleManager, 1, 0));
    
	NApiCore::ICustomPropertyManagerApi_1_0* particleCustomPropertyManager
		= static_cast<NApiCore::ICustomPropertyManagerApi_1_0*>(apiManager.getApi(eParticleCustomPropertyManager, 1, 0));

	iMAX_SCALE = particleCustomPropertyManager->getPropertyIndex(MAX_SCALE.c_str());
	
	// get the property index for the custom property Edit: Abhilash
	iCycle_Number = particleCustomPropertyManager->getPropertyIndex(Cycle_Number.c_str());

	return true;
}

ECalculateResult CParticleCycling::externalForce(
										int          threadId,
										double       time,
										double       timestep,
										int          id,
										const char   type[],
										double       mass,
										double       volume,
										double       density,
										unsigned int surfaces,
										double       posX,
										double       posY,
										double       posZ,
										double       velX,
										double       velY,
										double       velZ,
										double       angVelX,
										double       angVelY,
										double       angVelZ,
										const double orientation[9],
										NApiCore::ICustomPropertyDataApi_1_0* particlePropData,
										NApiCore::ICustomPropertyDataApi_1_0* simulationPropData,
										double&      calculatedForceX,
										double&      calculatedForceY,
										double&      calculatedForceZ,
										double&      calculatedTorqueX,
										double&      calculatedTorqueY,
										double&      calculatedTorqueZ)
{
	int rank;
	vector<string>::iterator it;

	it = find(m_particleNames.begin(), m_particleNames.end(), type);
	if (it != m_particleNames.end())
	{

		// New property of particle to store the cycle number during swell deswell :Edit : Abhi
		const double* CycleIndex = particlePropData->getValue(iCycle_Number);

		rank = (int)(it - m_particleNames.begin());

		// use dwell time to know the duration of swelling and deswelling of each particle Edit: Abhi
		double startTime = m_expansionStartTime[rank];
		double dwellTime = m_DwellTime[rank];

		if (time > startTime)
		{
			
			double ExpansionType	= m_expansionType[rank];
			double finalScale		= m_expansionFinalScale[rank];
			
			const double* particleMaxScale = particlePropData->getValue(iMAX_SCALE);
			if (particleMaxScale[0] == 0)
			{
				double* particleMaxScaleDelta = particlePropData->getDelta(iMAX_SCALE);
				particleMaxScaleDelta[0] = finalScale * m_particleMngr->getScale(id);
			};

			double currentScale = m_particleMngr->getScale(id);
			double scale = currentScale;

			// modifying the scaling wrt time step  : Edit :Abhilash
			// added an else statement to reduce the scale
			double cycle_position = fmod(time-startTime, (2 * dwellTime));

			if (cycle_position < dwellTime)
			{
				double fac = cycle_position / dwellTime; // Relative time for swelling

				if (ExpansionType == 1) //Expansion scale estimation for graphite
				{
					// delay the swelling of graphite untill 0.2 relative time
					// to reach 10% volume expansion at the end of swelling
					scale = 1 + (0.0375 * (fac-0.2));

					// Filter the values below 1 in the initial stages
					if (scale <= 1) 
					{
						scale = 1;
					}
				}
				
				if (ExpansionType == 2) //Expansion scale estimation for silicon
				{
					// start swelling of silicon from the begining
					// to reach 100% volume expansion at the end of swelling
					scale = 1 + (0.51125 * fac);
				}

				//stop swelling if any of the above cases reach maximum scale
				if (scale > particleMaxScale[0] && particleMaxScale[0] != 0)
				{
					scale = particleMaxScale[0];
				}
				
			}
			else
			{
				double fac = 1 - (cycle_position - dwellTime) / dwellTime; //	Relative time for de-swelling

				if (ExpansionType == 1) //Expansion scale estimation for graphite
				{
					// early start of graphite de-swelling
					scale = 1 + (0.0375 * (fac-0.2));
				}

				if (ExpansionType == 2) //Expansion scale estimation for silicon
				{
					scale = 1 + (0.51125 * fac);
					//do not de-swell unttill if the scale estimation is above the maximum scale
					if (scale > particleMaxScale[0] && particleMaxScale[0] != 0) 
					{
						scale = particleMaxScale[0];
					}
				}
				
				if (scale <= 1)
				{
					scale = 1;
				}
			}

			
			m_particleMngr->setScale(id, scale);
			
			// check the condition and increment the value of cycle index : Edit : Abhilash
			if ((cycle_position == 0.000) && (CycleIndex[0] <= m_maxcycles))
			{
				double* CycleIndexDelta = particlePropData->getDelta(iCycle_Number);
				CycleIndexDelta[0] = 1;
			}
		}

	}

    return eSuccess;
}

unsigned int CParticleCycling::getNumberOfRequiredProperties(
	const NApi::EPluginPropertyCategory category)
{
	if (eParticle == category)
	{
		return 2; // change the number of custom properties required Edit: Abhi
	}
	else
	{
		return 0;
	}
}

bool CParticleCycling::getDetailsForProperty(
	unsigned int                    propertyIndex,
	NApi::EPluginPropertyCategory   category,
	char                            name[NApi::CUSTOM_PROP_MAX_NAME_LENGTH],
	NApi::EPluginPropertyDataTypes& dataType,
	unsigned int&                   numberOfElements,
	NApi::EPluginPropertyUnitTypes& unitType,
	char                            initValBuff[NApi::BUFF_SIZE])
{
	if (0 == propertyIndex &&
		eParticle == category)
	{
		strcpy(name, MAX_SCALE.c_str());
		dataType = eDouble;
		numberOfElements = 1;
		unitType = eNone;

		std::ostringstream oss;
		oss << 0;
		strcpy(initValBuff, oss.str().c_str());

		//if you want to initialise the custom pperties with anything other than zero use below
		//example for csutom property with 3 elements

		//std::ostringstream oss;
		//oss << 1.0 << NApi::delim() << 0.0 << NApi::delim() << 1.0;
		//strcpy(initValBuff, oss.str().c_str());

		return true;
	}
	if (1 == propertyIndex &&
		eParticle == category)
	{
		strcpy(name, Cycle_Number.c_str());
		dataType = eDouble;
		numberOfElements = 1;
		unitType = eNone;

		std::ostringstream oss;
		oss << 0;
		strcpy(initValBuff, oss.str().c_str());
		return true;
	}
	else
	{
		return false;
	}
}