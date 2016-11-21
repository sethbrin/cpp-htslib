#include "easehts/args_parser.h"

namespace ncic {
namespace easehts {

std::vector<std::string>
ArgsParser::parse(int argc, char** argv) {

    std::vector<std::string> otherArguments;

    // first argument is the program name
    if ( argc >= 1 )
        programName = argv[0];

    // now, start iterating over each argument
    for(int argNumber=1; argNumber < argc; ++argNumber) {

        std::string argument = argv[argNumber];

        // arguments start with "-"
        // if not, push it into "other inputs"
        if ( argument.empty() )
            continue;

        if ( argument[0] != '-' ) {
            // add it as other argument and continue with the next arg
            otherArguments.push_back( argument );
            continue;
        }

        // this is a malformed argument:
        // "-"
        if ( argument.length() < 2 ) {
            std::stringstream error;
            error << "Malformed argument! (see arg number " << argNumber << ")";
            usage(error.str());
        }

        BaseOption* option = NULL;

        std::string possibleValue;

        // now, if the next char is a '-' it's a long option,
        // if not, it's a short one
        if ( argument[1] != '-' ) {

            option = findOption(argument[1]);

            // this looks like a short option, so let's check if there
            // are no more chars here, then we pick the value from here
            if ( argument.length() > 2 ) {
                possibleValue = argument.substr(2);
            }


        }
        else {
            // this looks like a long option, so let's check if there
            // are no more chars here, if not, that's malformed
            if ( argument.length() < 3 ) {
                std::stringstream error;
                error << "Malformed argument! (see arg number " << argNumber << ")";
                usage(error.str());
            }

            // let's allow the separation between key and value by '='
            // on long options
            std::string optionAndValueStr = argument.substr(2);

            // initially we suppose there is no value
            std::string optionStr = optionAndValueStr;

            int separator = optionStr.find('=');

            if ( separator != std::string::npos ) {
                optionStr     = optionAndValueStr.substr(0, separator);
                possibleValue = optionAndValueStr.substr(separator+1);
            }

            option = findOption( optionStr );

        }

        if ( option == NULL ) {
            std::stringstream error;
            error << "Unknown option '" << argument << "' (see arg number " << argNumber << ")";
            usage(error.str());
        }

        // let's see if this needs an argument
        if ( option->needArgument() ) {

            if ( possibleValue.empty() ) {

                // try to get the next one or fail
                if ( argNumber+1 < argc ) {

                    // let's move to the next argument
                    argNumber++;

                    option->setValue( argv[argNumber] );

                }
                else {
                    std::stringstream error;
                    error << "Option '" << argument << "' needs an additional argument";
                    usage(error.str());
                }

            }
            else {
                // the value we got directly from the option:
                //   --key=value or -kvalue
                option->setValue( possibleValue.c_str() );
            }

        }
        else {
            // set as read
            option->markAsFound();
        }

    }

    // now, let's do some basic checking

    // was the help option requested?
    if ( helpOption.isSet() ) {
        usage();
    }

    // let's go thru all the options to get all the
    // ones that are mandatories and that weren't
    // set
    std::string mandatoriesError;

    for(std::vector<BaseOption*>::iterator iter = options.begin();
        iter != options.end();
        ++iter
    ) {

        BaseOption* option = *iter;

        if ( option->isMandatory() && !option->isSet() ) {

            if ( ! mandatoriesError.empty() )
                mandatoriesError += ", ";

            mandatoriesError += getSummaryOptionText(option);
        }

    }

    if ( ! mandatoriesError.empty() ) {
        usage("The following arguments are mandatory: " + mandatoriesError);
    }

    return otherArguments;
}


void
ArgsParser::usage(const char* text) {

    if ( strcmp(text, "") != 0 ) {
        std::cerr << text << std::endl;
    }

    // add the options
    std::cerr << "Usage: ";
    std::cerr << programName << " ";

    std::stringstream optionsSummary;
    std::stringstream fullDescription;

    int maxWidth = 30;


    for(std::vector<BaseOption*>::iterator iter = options.begin();
        iter != options.end();
        ++iter
    ) {

        BaseOption* option = *iter;

        // this is the syntax:
        //   [ ] => optional
        //   short|long
        std::string summaryOptionBase = getSummaryOptionText(option);
        std::string fullOptionBase    = getFullOptionText(option);

        if ( option->needArgument() ) {
            summaryOptionBase += " value";
            fullOptionBase    += " value";
        }

        // summary
        if ( ! option->isMandatory() ) {
            optionsSummary << "[" << summaryOptionBase << "]";
        }
        else {
            optionsSummary << summaryOptionBase;
        }

        optionsSummary << " ";

        // full description
        fullDescription << " " << std::setw(maxWidth) << std::left << fullOptionBase << "\t\t" << option->getDescription() << std::endl;

    }

    std::cerr << optionsSummary.str() << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << fullDescription.str();

    // end
    exit(1);

}

std::string
ArgsParser::getOptionText(BaseOption* option, const char* separator) {

    std::string optionBase;

    if ( option->hasShortOption() ) {
        optionBase += "-";
        optionBase += option->getShortOption();
    }

    if ( option->hasLongOption() ) {

        if ( option->hasShortOption() )
            optionBase += separator;

        optionBase += "--" + option->getLongOption();
    }

    return optionBase;
}

BaseOption*
ArgsParser::findOption(char shortOption) {

    // iterate over the array and search for the short option
    for(int index=0; index < options.size(); ++index) {

        BaseOption* option = options.at(index);

        if ( option->matches( shortOption ) )
            return option;

    }

    return NULL;

}

BaseOption*
ArgsParser::findOption(std::string longOption) {

    // iterate over the array and search for the exact match option
    for(int index=0; index < options.size(); ++index) {
        BaseOption* option = options.at(index);
        if ( option->matches( longOption ) )
            return option;
    }

    // now, let's search for better matching ones...
    int bestMatchSize   = 0;
    BaseOption* bestMatchOption = NULL;
    std::vector<std::string> ambiguousOptions;

    for(int index=0; index < options.size(); ++index) {
        BaseOption* option = options.at(index);

        int matchSize = option->bestMatch(longOption);

        // we have a new winner
        if ( matchSize > bestMatchSize ) {
            bestMatchSize   = matchSize;
            bestMatchOption = option;

            // clear the list
            ambiguousOptions.clear();
        }
        else if (bestMatchSize>0 && matchSize==bestMatchSize ) {
            // this is conflicting with other option
            ambiguousOptions.push_back(option->getLongOption());
        }

    }

    if ( ambiguousOptions.size() > 0 ) {
        // if we're in an ambiguous case, it's better to report it
        std::stringstream error;
        error << "Option '" << longOption << "' is ambiguous: ";

        error << bestMatchOption->getLongOption();

        for(int index=0; index < ambiguousOptions.size(); ++index) {
            error << ", " << ambiguousOptions[index];
        }

        usage(error.str());
    }

    return bestMatchOption;

}

} // easehts
} // ncic
