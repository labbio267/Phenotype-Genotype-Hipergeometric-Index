# network_association
This code have been developed for the article: "Phenotype-loci associations in networks of patients with rare disorders: application to assist in the diagnosis of novel clinical cases. Anibal Bueno, Rocío Rodríguez-López, Armando Reyes-Palomares, Elena Rojano, Manuel Corpas, Julián Nevado, Pablo Lapunzina, Francisca Sánchez-Jiménez & Juan A.G. Ranea, European Journal of Human Genetics (2017), submitted."
This code has been created by Aníbal Bueno and modified by Elena Rojano.

The python code can be executed using the following structure:

python network_metrics.py -i --input_file -m --metric_type -o -output_file

Metric type can be one of the following:
	"hypergeometric"
	"jaccard"
	"PCC"
	"simpson"

Example:
python network_metrics.py -i example.txt -m hypergeometric -o example_output.txt

There is also provided an example network on which to apply the method, and the result after applying it (using Hypergeometric index):
example.txt
example_output.txt
