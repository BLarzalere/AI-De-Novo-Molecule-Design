# AI Molecule Design
A generative AI deep learning model for COVID-19 ligand design and screening using Python, Google Colab &amp; AutoDock Vina.

This repository contains the code and data used in the <a href="" target="_blank">AI Molecule Design</a> blog on Medium.com.  Provided below is a synopsis of the article.

The approach is to create an AI deep learning neural network that will learn how to create ligand molecules. A generative model is built from the trained neural network to design new synthetic molecules. AutoDock Vina is used to perform virtual screening of the molecules to assess their effectiveness for bonding to the COVID-19 main protease. The first generation of AI molecules is then used as the transfer learning dataset to further refine the neural network model for creating of more effective inhibitors of the COVID-19 main protease.

The high-level workflow is as follows:
- Obtain and prepare the molecule training data
- Create a neural network model and train it to make molecules
- Develop a generative AI model and create new molecules
- Virtually screen the new molecules against the COVID-19 main protease
- Use transfer learning to train the neural network on a more specific task
- Generate more AI molecules using the generative model
- Virtually screen the new molecules against the COVID-19 main protease
