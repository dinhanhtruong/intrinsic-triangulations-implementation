#include <QCoreApplication>
#include <QCommandLineParser>
#include <QDir>

#include <iostream>

#include "pathtracer.h"
#include "scene/scene.h"

#include <QImage>

#include "util/CS123Common.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addPositionalArgument("scene", "Scene file to be rendered");
    parser.addPositionalArgument("output", "Image file to write the rendered image to");
    parser.addPositionalArgument("renderVideo", "Render all frames for video output, true or false");
    parser.addPositionalArgument("algo", "Which algorithm to run if any: flip, refine");
    parser.addPositionalArgument("iterations", "Max number of flips/insertions for flip/refine respectively");
    parser.addPositionalArgument("minAngle", "Minimum angle bound for refinement in degrees");

    parser.process(a);

    const QStringList args = parser.positionalArguments();
    if(args.size() < 3) {
        std::cerr << "Error: Wrong number of arguments" << std::endl;
        a.exit(1);
        return 1;
    }
    QString scenefile = args[0];
    QString output = args[1];
    bool video = args[2] == "true" ? true : false;
    QString algo = args.size() > 3 ? args[3] : "";
    int iterations = args.size() > 4 ? args[4].toInt() : 0;
    double minAngle = args.size() > 5 ? args[5].toDouble() * M_PI / 180.0 : 0.0;

    QImage image(IMAGE_WIDTH, IMAGE_HEIGHT, QImage::Format_RGB32);

    Scene *scene;
    if(!Scene::load(scenefile, &scene)) {
        std::cerr << "Error parsing scene file " << scenefile.toStdString() << std::endl;
        a.exit(1);
        return 1;
    }

    PathTracer tracer(IMAGE_WIDTH, IMAGE_HEIGHT);

    QRgb *data = reinterpret_cast<QRgb *>(image.bits());

    if (video && (algo == "flip" || algo == "refine")) { // if no algo specified just render single still image of original mesh
        // if rendering a video need to get directory to save images
        QDir dir(output);
        if (!dir.exists()){
            // create directory if it doesn't exist
          dir.mkpath(".");
        }
        else {
            // clear any existing images in the directory
            dir.setNameFilters(QStringList() << "*.png");
            dir.setFilter(QDir::Files);
            for(QString& dirFile: dir.entryList()) {
                dir.remove(dirFile);
            }
        }
    }
    else {
        // for a single image add png extension to render single file
        if (!output.endsWith(".png")) output += ".png";
    }

    scene->getSPMesh()->setRenderInfo(scene, &tracer, &image, output, video);

    if (algo == "flip") {
        scene->getSPMesh()->renderFlipping(iterations);
    }
    else if (algo == "refine") {
        scene->getSPMesh()->renderRefine(minAngle, iterations);
    }
    else {
        // if no algo specified render single image
        scene->getSPMesh()->renderImage(output);
    }

    delete scene;
    a.exit();
}
