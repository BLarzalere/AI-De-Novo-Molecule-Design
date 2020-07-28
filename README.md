# AI--De-Novo-Molecule-Design
A generative AI deep learning model for COVID-19 ligand design and screening using Python, Google Colab &amp; AutoDock Vina

This repository contains the code and data used in the <a href="" target="_blank">AI - De Novo Molecule Design</a> blog on Medium.com.  Provided below is a synopis of the article.

The solution approach is to create an AI deep learning neural network that will learn how to create general molecules. Transfer learning is then used to teach the neural network to make molecules more specific to those needed to act as inhibitors to COVID-19 virus. A generative model is be built from the trained neural network to design new synthetic molecules. AutoDock Vina is used to perform virtual screening of the molecules to assess their effectiveness for bonding to the COVID-19 main protease.

The high-level workflow is as follows:
- Obtain and prepare the molecule training data
- Create a neural network model and train it to make molecules
- Use transfer learning to train the network on a more specific task
- Develop a generative AI model and create new molecules
- Virtually screen the new molecules against the COVID-19 main protease